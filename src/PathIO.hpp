#pragma once

// System-encoding adaptation for narrow-character file paths on Windows.
// The UI layer wxString::ToUTF8() gives UTF-8 bytes, while command-line argv / wxString::ToStdString()
// give ANSI (936/GBK) bytes, and the char* path of std::ifstream/ofstream is interpreted only according to the system ACP,
// so directories/file names containing Chinese cannot be opened. Here convert to UTF-16 with UTF-8 first and ACP fallback,
// then open via MSVC's wide-character stream overloads; pure ASCII paths take the original path, with unchanged behaviour.

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include <fstream>
#include <ios>
#include <string>

namespace pathio {

inline bool hasNonAscii(const std::string& s) {
    for (unsigned char c : s)
        if (c >= 0x80) return true;
    return false;
}

#ifdef _WIN32
inline std::wstring toWide(const std::string& path) {
    auto conv = [&path](UINT cp) -> std::wstring {
        const DWORD flags = (cp == CP_UTF8) ? MB_ERR_INVALID_CHARS : 0;
        const int n = MultiByteToWideChar(cp, flags, path.c_str(), -1, nullptr, 0);
        if (n <= 0) return {};
        std::wstring w(static_cast<size_t>(n - 1), L'\0');
        if (MultiByteToWideChar(cp, flags, path.c_str(), -1, w.data(), n) <= 0) return {};
        return w;
    };
    std::wstring w = conv(CP_UTF8);
    if (w.empty()) w = conv(CP_ACP);
    return w;
}
#endif

inline std::ifstream input(const std::string& path,
                           std::ios_base::openmode mode = std::ios::in) {
#ifdef _WIN32
    if (hasNonAscii(path)) {
        const std::wstring w = toWide(path);
        if (!w.empty()) {
            std::ifstream f(w.c_str(), mode);
            if (f) return f;
        }
    }
#endif
    return std::ifstream(path, mode);
}

inline std::ofstream output(const std::string& path,
                            std::ios_base::openmode mode = std::ios::out) {
#ifdef _WIN32
    if (hasNonAscii(path)) {
        const std::wstring w = toWide(path);
        if (!w.empty()) return std::ofstream(w.c_str(), mode);
    }
#endif
    return std::ofstream(path, mode);
}

}  // namespace pathio
