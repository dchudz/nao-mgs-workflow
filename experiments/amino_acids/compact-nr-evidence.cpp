// Preserve the sufficient statistic used by protein-nr-summary.py.
// Raw alignments stay untouched. Control and seed subsets remain raw rows.
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;
struct Best { double score; std::string score_text; };
uint64_t number(std::string_view s) {
    uint64_t n = 0;
    auto r = std::from_chars(s.data(), s.data() + s.size(), n);
    if (r.ec != std::errc() || r.ptr != s.data() + s.size())
        throw std::runtime_error("Invalid unsigned integer field");
    return n;
}
double score_value(std::string_view s) {
    double n = 0;
    auto r = std::from_chars(s.data(), s.data() + s.size(), n);
    if (r.ec != std::errc() || r.ptr != s.data() + s.size() || !std::isfinite(n))
        throw std::runtime_error("Invalid finite score field");
    return n;
}
bool digits(std::string_view s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; });
}
int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("Usage: compact-nr-evidence INPUT NEW_OUTPUT_DIRECTORY EXPECTED_ROWS_OR_0");
    fs::path input(argv[1]), out(argv[2]);
    uint64_t expected = number(argv[3]);
    if (!fs::create_directory(out)) throw std::runtime_error("Output directory must be new");
    std::vector<char> buffer(8 << 20);
    std::ifstream in;
    in.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
    in.open(input, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open input");
    std::ofstream compact(out / "query-calls.tsv.partial", std::ios::binary);
    std::ofstream controls(out / "native-control-alignments.tsv.partial", std::ios::binary);
    std::ofstream seeds(out / "seed-alignments.tsv.partial", std::ios::binary);
    for (auto* f : {&compact, &controls, &seeds}) f->exceptions(std::ios::badbit | std::ios::failbit);
    uint64_t rows = 0, bytes = 0, eligible = 0, emitted = 0, control_rows = 0, seed_rows = 0;
    std::string line, query;
    std::unordered_set<std::string> seen;
    std::unordered_map<uint64_t, Best> best;
    bool control = false, seed = false;
    auto flush = [&]() {
        if (best.empty()) return;
        auto maximum = std::max_element(best.begin(), best.end(), [](const auto& a, const auto& b) { return a.second.score < b.second.score; });
        std::vector<uint64_t> taxa;
        for (const auto& x : best) if (x.second.score >= 0.9 * maximum->second.score) taxa.push_back(x.first);
        std::sort(taxa.begin(), taxa.end());
        compact << query << '\t' << maximum->second.score_text << '\t';
        for (size_t i = 0; i < taxa.size(); ++i) { if (i) compact << ';'; compact << taxa[i]; }
        compact << '\n'; ++emitted;
        best.clear();
    };
    while (std::getline(in, line)) {
        if (in.eof()) throw std::runtime_error("Input must end with a complete newline");
        ++rows; bytes += line.size() + 1;
        std::array<std::string_view, 7> f;
        size_t start = 0;
        for (size_t i = 0; i < 6; ++i) {
            size_t tab = line.find('\t', start);
            if (tab == std::string::npos) throw std::runtime_error("Expected seven fields");
            f[i] = std::string_view(line).substr(start, tab - start); start = tab + 1;
        }
        f[6] = std::string_view(line).substr(start);
        if (f[6].find('\t') != std::string_view::npos || f[6].find('\r') != std::string_view::npos)
            throw std::runtime_error("Unexpected delimiter or carriage return");
        if (f[0] != query) {
            flush(); query.assign(f[0]);
            if (query.empty() || !seen.insert(query).second)
                throw std::runtime_error("Query groups must be contiguous and nonempty");
            control = query.rfind("NR_CONTROL_", 0) == 0;
            seed = query.rfind("SEED_", 0) == 0;
        }
        if (control) { controls << line << '\n'; ++control_rows; }
        if (seed) { seeds << line << '\n'; ++seed_rows; }
        if (number(f[4]) < 25) continue;
        ++eligible;
        double score = score_value(f[2]);
        auto retain = [&](uint64_t taxon) {
            auto old = best.find(taxon);
            if (old == best.end() || score > old->second.score)
                best.insert_or_assign(taxon, Best{score, std::string(f[2])});
        };
        bool any = false;
        size_t pos = 0;
        while (true) {
            size_t sep = f[6].find(';', pos);
            auto token = f[6].substr(pos, sep == std::string_view::npos ? sep : sep - pos);
            if (digits(token)) { retain(number(token)); any = true; }
            if (sep == std::string_view::npos) break;
            pos = sep + 1;
        }
        if (!any) retain(0);
    }
    if (!in.eof()) throw std::runtime_error("Input read failed");
    flush();
    if (bytes != fs::file_size(input) || (expected && rows != expected))
        throw std::runtime_error("Input byte or expected row count mismatch");
    for (auto* f : {&compact, &controls, &seeds}) f->close();
    for (const auto& name : {"query-calls.tsv", "native-control-alignments.tsv", "seed-alignments.tsv"})
        fs::rename(out / (std::string(name) + ".partial"), out / name);
    std::ofstream stats(out / "counts.json");
    stats.exceptions(std::ios::badbit | std::ios::failbit);
    stats << "{\"input_rows\":" << rows << ",\"input_bytes\":" << bytes
          << ",\"aligned_query_mates\":" << seen.size() << ",\"eligible_rows\":" << eligible
          << ",\"near_best_query_rows\":" << emitted << ",\"native_control_raw_rows\":" << control_rows
          << ",\"seed_raw_rows\":" << seed_rows << "}\n";
    stats.close();
    std::cout << "COMPLETE input_rows=" << rows << " compact_rows=" << emitted << '\n';
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}
