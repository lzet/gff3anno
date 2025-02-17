#ifndef GETOPTS_H
#define GETOPTS_H
#include <string>
#include <vector>

/**
 * @brief The cmd_opt_t class
 * поле name - входной парамерт без начального символа '-'
 * поле values - параметры, стоящие за параметром name
 */
struct cmd_opt_t {
    std::string name;
    std::vector<std::string> values;
    /**
     * @brief equal соответствует ли имя параметра заданному имени optname
     * @param optname имя, с которым сравниваем
     * @param canstartpart имя может начинаться с optname и в конце содержать любые символы
     * @return true если совпадает
     */
    bool equal(const std::string &optname, bool canstartpart = false) const;
};

class get_opts_t
{
    std::string _error;
    std::string _program_name;
    std::string _program_path;
    std::vector<cmd_opt_t> _opts;
public:
    get_opts_t(int argc, char **argv);
    /**
     * @brief result
     * @return результат разбора коммандной строки
     */
    std::vector<cmd_opt_t>& result();
    /**
     * @brief program_name
     * @return имя программы взятое из argv[0]
     */
    std::string& program_name();
    /**
     * @brief program_path
     * @return абсолютный путь до программы из argv[0]
     */
    std::string& program_path();
    /**
     * @brief dump
     * @param prefix префикс для каждой строки
     * @return result в виде форматированного текста (для вывода в std::cout, к примеру)
     */
    std::string dump(const std::string &prefix = "") const;
};

#endif // GETOPTS_H
