/* MIT License
 *
 * Copyright (c) 2020 Aleksa Ilic <aleksa.d.ilic@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#include <algorithm>
#include <iostream>
#include <numeric>
#include <regex>
#include <string>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/optional/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/rational.hpp>

#ifdef TESTING_ENABLED
#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <prettyprint.hpp>
#else
#endif

#include <fort.hpp>

static constexpr auto kVersion = "v0.1.0";
static constexpr auto kProgramName = "OPNA_2: Sturm's Theorem Evaluator";
static constexpr auto kUnderlineType = '-';

// -- CHANGE THESE TYPES FOR MORE PRECISION
namespace mp = boost::multiprecision;
using IntType = mp::int1024_t;
using FloatType = mp::number<mp::cpp_dec_float<256>>;
using Polyonimal = boost::math::tools::polynomial<IntType>;

enum class Field {

};

constexpr const char *FieldToString(Field field) {
    switch (field) {
        default:
            return "Unknown";
    }
}

struct Config {
    size_t iterations = 14;
    const struct ft_border_style *table_style = FT_NICE_STYLE;
    std::vector<Field> displayed_fields = {};

};
static const Config kDefaultConfig;

/// Performs substring operation with range check without throwing errors
/// \param string String to be processed
/// \param position Position from which to extract
/// \param length How many characters are to be processed
/// \return Extracted substring if range check passes or boost::none
boost::optional<std::string>
SafeSubstring(const std::string &string, size_t position, size_t length = std::string::npos) {
    if (position >= string.length())
        return boost::none;
    else
        return string.substr(position, length);
}

/// Retrieves program's base name from the exec path
static const char *GetProgramName(const char *path) {
    const char *last = path;
    while (*path++) {
        last = (*path == '/' || *path == '\\') ? path : last;
    }
    return last;
}

#ifdef TESTING_ENABLED
TEST_CASE("GetProgramName") {
    constexpr auto kExamplePath1 = "/home/ilic/opna-1/opna_1";
    constexpr auto kExamplePath2 = "C:\\Users\\My Documents\\opna_1";
    constexpr auto kExpectedName1 = "/opna_1";
    constexpr auto kExpectedName2 = "\\opna_1";

    REQUIRE(strcmp(GetProgramName(kExamplePath1), kExpectedName1) == 0);
    REQUIRE(strcmp(GetProgramName(kExamplePath2), kExpectedName2) == 0);
}
#else

int main(int argc, const char *argv[]) {
    namespace po = boost::program_options;

    Config config = kDefaultConfig;
    auto program_name = GetProgramName(argv[0]);
    std::string polynomial;
    std::string table_style;
    std::string fields;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "print help message")
            ("version", "print version information")
            ("examples", "show examples")
            ("table", "print evaluation table of every iteration")
            ("table_style", po::value<std::string>(&table_style), "Specify table style: nice|double|simple|empty")
            ("fields", po::value<std::string>(&fields),
             "Comma delimited list of fields to be shown: ")
            ("polynomial", po::value<std::string>(&polynomial), "polynomial to be parsed");

    po::positional_options_description pos_desc;
    pos_desc.add("polynomial", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "Usage: ." << program_name << " [OPTIONS] <POLYNOMIAL> \n"
                  << desc << '\n';
        return 0;
    }

    if (vm.count("version")) {
        std::ostringstream name_header;
        name_header << kProgramName << ' ' << kVersion;

        std::string underline(name_header.str().size(), kUnderlineType);

        std::cout << name_header.str() << "\n"
                  << underline << "\n"
                  << "Integer type w/: " << std::numeric_limits<IntType>::digits << "bits\n";
        return 0;
    }

    if (vm.count("examples")) {
        fort::char_table table;
        table.set_border_style(FT_EMPTY_STYLE);

        const std::vector<std::tuple<std::string, std::string>> examples = {
                std::make_tuple("example_name", "example_description")
        };

        table << fort::header
              << "" << "Command" << "Description" << fort::endr;

        for (size_t i = 0; i < examples.size(); i++) {
            const auto&[command, description] = examples[i];
            table << i + 1 << std::string(".") + program_name + " " + command << description << fort::endr;
        }

        std::cout << table.to_string();
        return 0;
    }

    if (!vm.count("polynomial")) {
        std::cerr << "You need to specify a polynomial. Use --help for usage information \n";
        return 1;
    }

    if (vm.count("table_style")) {
        if (boost::iequals(table_style, "nice")) {
            config.table_style = FT_NICE_STYLE;
        } else if (boost::iequals(table_style, "simple")) {
            config.table_style = FT_SIMPLE_STYLE;
        } else if (boost::iequals(table_style, "double")) {
            config.table_style = FT_DOUBLE2_STYLE;
        } else if (boost::iequals(table_style, "empty")) {
            config.table_style = FT_EMPTY_STYLE;
        } else {
            std::cerr << "Incorrect table style supplied: " << table_style << ". Use --help for usage information \n";
            return 1;
        }
    }

    if (vm.count("fields")) {
        config.displayed_fields.clear();
        std::vector<std::string> field_vector;
        boost::split(field_vector, fields, [](char c) { return c == ','; });
        for (const auto &field:field_vector) {
            std::cerr << "Incorrect field name supplied: " << field << ". Use --help for usage information \n";
            return 1;
        }
    }
}

#endif