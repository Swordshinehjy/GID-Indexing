// ============================================================================
//  parallel.h — minimal parallel for: split [0,n) into chunks handed to std::thread.
//  Used only in dependency-free candidate enumeration ((u,v) pages, l-triple shards, random-trial shards).
// ============================================================================
#pragma once

#include <algorithm>
#include <atomic>
#include <functional>
#include <thread>
#include <vector>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#endif

namespace gid {

// Physical core count. std::thread::hardware_concurrency() gives the **logical processor** count
// (including hyperthreading / big.LITTLE SMT siblings); spawning threads by it makes two threads
// contend for one physical core, which is measurably slower (~10-30% for memory/cache-bound kernels). Falls back to the logical count on failure.
inline unsigned int physicalCores() {
#if defined(_WIN32)
    // Each RelationProcessorCore record = one physical core (possibly with several logical processor masks)
    DWORD len = 0;
    if (!GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &len) &&
        GetLastError() != ERROR_INSUFFICIENT_BUFFER)
        return std::thread::hardware_concurrency();
    std::vector<char> buf(len);
    if (!GetLogicalProcessorInformationEx(
            RelationProcessorCore,
            reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf.data()), &len))
        return std::thread::hardware_concurrency();
    unsigned int count = 0;
    for (DWORD off = 0; off < len;) {
        const auto* e = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(
            buf.data() + off);
        if (e->Relationship == RelationProcessorCore) ++count;
        if (e->Size == 0) break;
        off += e->Size;
    }
    return count > 0 ? count : std::thread::hardware_concurrency();
#elif defined(__APPLE__)
    int n = 0;
    size_t sz = sizeof(n);
    if (sysctlbyname("hw.physicalcpu", &n, &sz, nullptr, 0) == 0 && n > 0)
        return (unsigned int)n;
    return std::thread::hardware_concurrency();
#elif defined(__linux__)
    // "cpu cores" in /proc/cpuinfo is the **number of cores per physical CPU**, so multiply by the
    // physical package count; the product is inaccurate for hybrid architectures but closer than the logical count; fall back on failure.
    unsigned int cores = 0, pkgs = 0;
    if (FILE* f = std::fopen("/proc/cpuinfo", "r")) {
        char line[256];
        while (std::fgets(line, sizeof(line), f)) {
            if (std::strncmp(line, "cpu cores", 9) == 0) {
                const char* c = std::strchr(line, ':');
                if (c && cores == 0) cores = (unsigned int)std::atoi(c + 1);
            } else if (std::strncmp(line, "physical id", 11) == 0) {
                const char* c = std::strchr(line, ':');
                if (c) pkgs = std::max(pkgs, (unsigned int)std::atoi(c + 1) + 1u);
            }
        }
        std::fclose(f);
    }
    if (cores > 0 && pkgs > 0) return cores * pkgs;
    return std::thread::hardware_concurrency();
#else
    return std::thread::hardware_concurrency();
#endif
}

// Default parallelism: physical core count (capped at 16 so thread-creation overhead does not
// dominate small tasks). The function-local static queries platform info once (thread-safe since C++11).
inline unsigned int hardwareThreads() {
    static const unsigned int n = [] {
        unsigned int c = physicalCores();
        if (c == 0) c = std::thread::hardware_concurrency();
        if (c == 0) c = 4;
        return std::min(c, 16u);
    }();
    return n;
}

// Dynamic task queue: task granularity 1, suitable when task durations vary widely (e.g.
// candidate-level Part TWO). Static chunking in parallelFor has a pronounced long tail when durations are uneven.
inline void parallelForEach(int n, const std::function<void(int)>& body, int threads = -1) {
    if (n <= 0) return;
    if (threads < 0) threads = (int)hardwareThreads();
    threads = std::min(threads, n);
    if (threads <= 1) {
        for (int i = 0; i < n; ++i) body(i);
        return;
    }
    std::atomic<int> next{0};
    std::vector<std::thread> ts;
    ts.reserve((size_t)threads);
    for (int t = 0; t < threads; ++t)
        ts.emplace_back([&]() {
            for (;;) {
                const int i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= n) break;
                body(i);
            }
        });
    for (auto& th : ts) th.join();
}

// body(begin,end); called serially when threads<=1
inline void parallelFor(int begin, int end, const std::function<void(int, int)>& body,
                        int threads = -1) {
    if (end <= begin) return;
    if (threads < 0) threads = (int)hardwareThreads();
    if (threads <= 1 || end - begin < 64) {
        body(begin, end);
        return;
    }
    const int total = end - begin;
    const int chunk = std::max(1, (total + threads - 1) / threads);
    std::vector<std::thread> ts;
    ts.reserve(threads);
    for (int s = 0; s < threads; ++s) {
        const int b = begin + s * chunk;
        const int e = std::min(end, b + chunk);
        if (b >= e) break;
        ts.emplace_back(body, b, e);
    }
    for (auto& t : ts) t.join();
}

} // namespace gid
