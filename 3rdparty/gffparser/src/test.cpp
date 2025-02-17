#include <gffparser.h>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace gffparser;

int main(int argc, char **argv) {
    if(argc < 2) {
        std::cerr << "USAGE: program path/to/file.gff3" << std::endl;
        return 1;
    }
    std::ifstream ifs(argv[1]);
    if(!ifs.is_open()) {
        std::cerr << "can't open file: " << argv[1] << std::endl;
        return 1;
    }
    std::string line;
    gff_parser_t gff(4);
    using std::chrono::steady_clock;

    steady_clock::time_point init_start = steady_clock::now();
    while(std::getline(ifs, line)) {
        gff << line;
        if(gff.has_error()) {
            std::cerr << gff.error() << std::endl;
            return 1;
        }
    }
    gff.flush();
    std::cout << "init time in sec: "
              << std::chrono::duration_cast<std::chrono::seconds>(steady_clock::now()-init_start).count()
              << std::endl;
    std::cout << "db size: " << gff.size() << std::endl;

    std::cout << "select type=gene gene_id=ENSG00000310526.1" << std::endl;
    init_start = steady_clock::now();
    auto gene_id = gff.get_by("gene", {gff_attribute_t("gene_id", "ENSG00000310526.1")});
    std::cout << "select by gene_id time in msec: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now()-init_start).count()
              << std::endl;
    std::cout << "RESULT1:" << std::endl;
    for(const auto &g: gene_id)
        std::cout << g->type << ": " << g->get_attr("gene_name").get_string() << " / " << g->get_attr("gene_id").get_string() << std::endl;
        //std::cout << g->str() << std::endl;

    std::cout << "select type=gene position=chr20:44619522-44651699" << std::endl;
    init_start = steady_clock::now();
    auto gene_pos = gff.get_by("gene", gff_postition_t("chr20", 44619522, 44651699));
    std::cout << "select by gene position time in msec: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now()-init_start).count()
              << std::endl;
    std::cout << "RESULT2:" << std::endl;
    for(const auto &g: gene_pos)
        std::cout << g->type << ": " << g->get_attr("gene_name").get_string() << " / " << g->get_attr("gene_id").get_string() << std::endl;
        //std::cout << g->str() << std::endl;

    return 0;
}
