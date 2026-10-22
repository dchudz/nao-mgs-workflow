// Copy every raw alignment for specified query IDs, preserving bytes and order.
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
namespace fs = std::filesystem;
int main(int argc, char** argv) {
  try {
    if (argc != 5) throw std::runtime_error("Usage: extract-nr-query-subset INPUT IDS NEW_DIRECTORY EXPECTED_ROWS_OR_0");
    fs::path input(argv[1]), ids(argv[2]), out(argv[3]);
    uint64_t expected = std::stoull(argv[4]), rows = 0, bytes = 0, selected = 0;
    std::unordered_set<std::string> wanted, found;
    std::string line;
    std::ifstream names(ids);
    if (!names) throw std::runtime_error("Cannot read query IDs");
    while (std::getline(names, line)) {
        if (line.empty() || line.find_first_of("\t\r ") != std::string::npos || !wanted.insert(line).second)
            throw std::runtime_error("Query IDs must be nonempty, distinct single tokens");
    }
    if (!fs::create_directory(out)) throw std::runtime_error("Output directory must be new");
    std::vector<char> buffer(8 << 20);
    std::ifstream in;
    in.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
    in.open(input, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read raw alignments");
    std::ofstream dest(out / "alignments.tsv.partial", std::ios::binary);
    dest.exceptions(std::ios::badbit | std::ios::failbit);
    std::string previous;
    bool keep = false;
    while (std::getline(in, line)) {
        if (in.eof()) throw std::runtime_error("Incomplete final input line");
        ++rows; bytes += line.size() + 1;
        auto tab = line.find('\t');
        if (tab == std::string::npos) throw std::runtime_error("Missing query delimiter");
        if (line.compare(0, tab, previous) != 0) {
            previous.assign(line, 0, tab); keep = wanted.count(previous);
            if (keep) found.insert(previous);
        }
        if (keep) { dest << line << '\n'; ++selected; }
    }
    if (!in.eof() || bytes != fs::file_size(input) || (expected && rows != expected))
        throw std::runtime_error("Input read/count mismatch");
    dest.close();
    fs::rename(out / "alignments.tsv.partial", out / "alignments.tsv");
    std::ofstream stats(out / "counts.json");
    stats.exceptions(std::ios::badbit | std::ios::failbit);
    stats << "{\"input_rows\":" << rows << ",\"input_bytes\":" << bytes
          << ",\"selected_raw_rows\":" << selected << ",\"requested_queries\":" << wanted.size()
          << ",\"aligned_selected_queries\":" << found.size() << "}\n";
    stats.close();
    std::cout << "EXTRACTED " << selected << " raw rows\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
