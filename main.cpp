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
#include <map>
#include <numeric>
#include <regex>
#include <string>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/optional.hpp>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/program_options.hpp>

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
using Fraction = boost::rational<IntType>;
using CoefficientType = Fraction;
using PowerType = unsigned long long;
using Polynomial = boost::math::tools::polynomial<CoefficientType>;

struct Config {
    const struct ft_border_style *table_style = FT_SIMPLE_STYLE;
};
static const Config kDefaultConfig;

/// Retrieves program's base name from the exec path
static const char *GetProgramName(const char *path) {
    const char *last = path;
    while (*path++) {
        last = (*path == '/' || *path == '\\') ? path : last;
    }
    return last;
}

std::string PolynomialToHumanString(const Polynomial &polynomial) {
    if (polynomial.size() == 0) {
        return std::string("");
    }
    std::ostringstream oss;
    auto max_power = polynomial.size() - 1;
    for (PowerType power = max_power; power >= 0 && power != std::numeric_limits<PowerType>::max(); power--) {
        const auto &coefficient = polynomial[power];
        if (coefficient != 0) {
            if (coefficient > 0 && power != max_power) {
                oss << '+';
            } else if (coefficient < 0) {
                oss << '-';
            }
            if (!(mp::abs(coefficient.numerator()) == 1 && mp::abs(coefficient.denominator()) == 1) || power == 0) {
                oss << mp::abs(coefficient.numerator());
                if (mp::abs(coefficient.denominator()) != 1) {
                    oss << '/' << mp::abs(coefficient.denominator());
                }
            }
            if (power > 1) {
                oss << "x^" << std::to_string(power);
            } else if (power == 1) {
                oss << "x";
            }
        }
    }
    return oss.str();
}

Polynomial GetDerivative(const Polynomial &polynomial) {
    if (polynomial.size() <= 1) {
        return Polynomial{};
    }

    auto derivative = Polynomial(std::next(polynomial.data().begin()), polynomial.data().end());
    for (PowerType power = 1; power <= derivative.size(); power++) {
        derivative[power - 1] = power * derivative[power - 1];
    }

    return derivative;
}

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

/// Parse decimal number string as fraction
/// \param number Number to be parsed
/// \return Corresponding fraction
Fraction DecimalToFraction(std::string number) {
    auto dot_index = number.find('.');
    if (dot_index == std::string::npos) {
        return Fraction(IntType(number), 1);
    }
    std::string before_point = number.substr(0, dot_index);
    std::string after_point = SafeSubstring(number, dot_index + 1).get_value_or("");

    auto numerator = IntType(before_point + after_point);
    auto denominator = mp::pow(IntType(10), after_point.size());

    return Fraction(numerator, denominator);
}

Fraction FractionFromString(const std::string &fraction_string) {
    std::smatch match;
    if (std::regex_match(fraction_string, match, std::regex("^([0-9.]+)/([0-9.]+)$"))) {
        assert(match.size() == 3);
        Fraction numerator = DecimalToFraction(match[1].str());
        Fraction denominator = DecimalToFraction(match[2].str());
        return numerator / denominator;
    } else if (std::regex_match(fraction_string, match, std::regex("^([0-9.]+)$"))) {
        assert(match.size() == 2);
        return DecimalToFraction(match[1].str());
    } else {
        throw std::invalid_argument("Ill-formatted fraction: " + fraction_string);
    }
}

Polynomial ParsePolynomial(std::string polynomial_string) {
    static const auto kVariableExtractionRegex = std::regex(
            "([+-]?(?:(?:[0-9/.]+x\\^[0-9/.]+)|(?:[0-9/.]+x)|(?:[0-9/.]+)|(?:x)))");
    static const auto kVariableWithoutPowerRegex = std::regex("^[+-]?([0-9/.]+)x$");
    static const auto kVariableWithPowerRegex = std::regex("^[+-]?([0-9/.]+)x\\^([0-9]+)$");
    static const auto kConstantRegex = std::regex("^[+-]?([0-9/.]+)$");

    for (size_t it = polynomial_string.find('x'); it != std::string::npos; it = polynomial_string.find('x', ++it)) {
        if (it == 0 || !std::isdigit(polynomial_string[it - 1])) {
            polynomial_string.insert(it++, "1");
        }
    }

    std::vector<std::string> variables;
    std::copy(std::sregex_token_iterator(polynomial_string.begin(), polynomial_string.end(), kVariableExtractionRegex),
              std::sregex_token_iterator(),
              std::back_inserter(variables));

    std::map<PowerType, CoefficientType> polynomial_map;
    for (const auto &variable : variables) {
        int sign = (variable.front() == '-' ? -1 : 1);

        std::smatch match;
        PowerType power;
        if (std::regex_match(variable, match, kVariableWithPowerRegex)) {
            assert(match.size() == 3);
            power = std::stoull(match[2].str());
        } else if (std::regex_match(variable, match, kVariableWithoutPowerRegex)) {
            assert(match.size() == 2);
            power = 1;
        } else if (std::regex_match(variable, match, kConstantRegex)) {
            assert(match.size() == 2);
            power = 0;
        } else {
            throw std::invalid_argument("Ill-formatted polynomial" + polynomial_string);
        }
        auto coefficient = FractionFromString(match[1].str());

        polynomial_map[power] += sign * coefficient;
    }

    if (!polynomial_map.empty()) {
        auto max_power = polynomial_map.rbegin()->first;
        std::vector<CoefficientType> polynomial_vector(max_power + 1);

        for (const auto&[power, coefficient] : polynomial_map) {
            polynomial_vector[power] = coefficient;
        }

        return Polynomial{polynomial_vector};
    }

    return Polynomial{};
}

std::vector<Polynomial> SturmSequence(const Polynomial &polynomial) {
    std::vector<Polynomial> iterations;
    iterations.push_back(polynomial);
    iterations.push_back(GetDerivative(polynomial));

    while (!iterations.back().is_zero()) {
        const auto &p0 = *std::next(iterations.rbegin());
        const auto &p1 = *iterations.rbegin();
        iterations.push_back(-1 * (p0 % p1));
    }
    iterations.pop_back();

    return iterations;
}

std::vector<bool> SignSequence(const std::vector<Polynomial> &sturm_sequence, IntType alpha) {
    std::vector<bool> sign_sequence;
    sign_sequence.reserve(sturm_sequence.size());
    std::transform(sturm_sequence.begin(), sturm_sequence.end(), std::back_inserter(sign_sequence),
                   [alpha](const auto &polynomial) {
                       return polynomial.evaluate(alpha) >= 0;
                   });
    return sign_sequence;
}

IntType SignVariations(const std::vector<bool> &sign_sequence) {
    if (sign_sequence.size() <= 1) {
        return sign_sequence.size();
    }
    IntType sign_variations = 0;
    for (auto sign_it = std::next(sign_sequence.begin()); sign_it != sign_sequence.end(); ++sign_it) {
        if (*std::prev(sign_it) != *sign_it) {
            sign_variations++;
        }
    }
    return sign_variations;
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

TEST_CASE("Polynomial conversions") {
    const auto kPolynomialTests = std::vector<std::pair<std::string, Polynomial>>{
            std::make_pair("x^4+x^3-x-1", Polynomial{-1, -1, 0, 1, 1}),
            std::make_pair("-x^4+x^3-x-1", Polynomial{-1, -1, 0, 1, -1}),
            std::make_pair("-1/2x^2", Polynomial{ 0/1, 0/1, Fraction(-1,2) }),
            std::make_pair("4x^3+3x^2-1", Polynomial{ -1, 0, 3, 4}),
            std::make_pair("3/16x^2+3/4x+15/16", Polynomial{ Fraction(15,16), Fraction(3,4), Fraction(3,16)}),
            std::make_pair("-32x-64", Polynomial{ -64, -32}),
            std::make_pair("-3/16", Polynomial{Fraction(-3,16)})
    };
    for (const auto &[polynomial_string, expected_result] : kPolynomialTests) {
        const auto parsed_polynomial = ParsePolynomial(polynomial_string);
        REQUIRE(expected_result == parsed_polynomial);
        REQUIRE(polynomial_string == PolynomialToHumanString(parsed_polynomial));
    }
}
TEST_CASE("Sturm iterations"){
    constexpr auto kPolynomialString = "x^4+x^3-x-1";
    constexpr auto kRangeFrom = -5;
    constexpr auto kRangeTo = 5;
    const auto kExpectedIterations = std::vector<Polynomial>{
        ParsePolynomial(kPolynomialString),
        ParsePolynomial("4x^3+3x^2-1"),
        ParsePolynomial("3/16x^2+3/4x+15/16"),
        ParsePolynomial("-32x-64"),
        ParsePolynomial("-3/16")
    };
    constexpr auto kExpectedEvaluationFrom = 3;
    constexpr auto kExpectedEvaluationTo = 1;
    auto iterations = SturmSequence(ParsePolynomial(kPolynomialString));
    REQUIRE(kExpectedIterations == iterations);
    REQUIRE(kExpectedEvaluationFrom == SignVariations(SignSequence(iterations, kRangeFrom)));
    REQUIRE(kExpectedEvaluationTo == SignVariations(SignSequence(iterations, kRangeTo)));
}

#else

std::string SignSequenceToString(const std::vector<bool> &sign_sequence) {
    std::ostringstream oss;
    oss << '(';
    for (const bool is_plus : sign_sequence) {
        oss << (is_plus ? '+' : '-') << ',';
    }
    if (!sign_sequence.empty()) {
        oss.seekp(-1, std::ios_base::end);
    }
    oss << ')';
    return oss.str();
}

int main(int argc, const char *argv[]) {
    namespace po = boost::program_options;

    Config config = kDefaultConfig;
    auto program_name = GetProgramName(argv[0]);
    std::string polynomial;
    std::string table_style;
    std::string fields;
    IntType range[2];

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "print help message")
            ("version", "print version information")
            ("examples", "show examples")
            ("verbose", "detailed output of every step")
            ("table_style", po::value<std::string>(&table_style),
             "Specify table style for verbose output: nice|double|simple|empty")
            ("polynomial", po::value<std::string>(&polynomial), "polynomial to be parsed")
            ("from", po::value<IntType>(&range[0]), "interval from")
            ("to", po::value<IntType>(&range[1]), "interval to");

    po::positional_options_description pos_desc;
    pos_desc.add("polynomial", 1);
    pos_desc.add("from", 1);
    pos_desc.add("to", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).positional(pos_desc).options(desc).style(
            po::command_line_style::unix_style ^ po::command_line_style::allow_short
    ).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "Usage: ." << program_name << " [OPTIONS] <POLYNOMIAL> <INTERVAL_FROM> <INTERVAL_TO> \n"
                  << desc << '\n';
        return 0;
    }

    if (vm.count("version")) {
        std::ostringstream name_header;
        name_header << kProgramName << ' ' << kVersion;

        std::string underline(name_header.str().size(), kUnderlineType);

        std::cout << name_header.str() << "\n"
                  << underline << "\n"
                  << "Coefficient = Integer w/: " << std::numeric_limits<IntType>::digits << "bits\n"
                  << "Power = Integer w/: " << std::numeric_limits<PowerType>::digits << "bits\n";
        return 0;
    }

    if (vm.count("examples")) {
        fort::char_table table;
        table.set_border_style(FT_EMPTY_STYLE);

        const std::vector<std::tuple<std::string, std::string>> examples = {
                std::make_tuple("x^4+x^3-x-1 -5 5", "Evaluate P(x) = x^4+x^3-x-1 with Sturm's theorem on range (-5,5]"),
                std::make_tuple("x^4+x^3-x-1 -5 5 --verbose", "Same as above but with detailed output of every step"),
                std::make_tuple("1/2x^2-5", "Fractions are supported as well"),
                std::make_tuple("x^3+x^5+2x^3", "Mixed ordering and repeating of the same element is allowed")
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

    const auto sturm_sequence = SturmSequence(ParsePolynomial(polynomial));
    const auto sign_sequences = std::vector<std::vector<bool>>{SignSequence(sturm_sequence, range[0]),
                                                               SignSequence(sturm_sequence, range[1])};
    const auto sign_variations = std::vector<IntType>{SignVariations(sign_sequences[0]),
                                                      SignVariations(sign_sequences[1])};
    const auto number_of_roots = sign_variations[0] - sign_variations[1];

    const auto print_sturm_sequence = [&]() {
        fort::char_table table;
        table.set_border_style(config.table_style);
        table << fort::header << "Sturm sequence" << fort::endr;
        for (size_t i = 0; i < sturm_sequence.size(); i++) {
            std::ostringstream oss;
            oss << "p" << i << "(x) = " << PolynomialToHumanString(sturm_sequence[i]);
            table << oss.str() << fort::endr;
        }
        std::cout << table.c_str();
    };

    const auto print_sign_table = [&]() {
        fort::char_table table;
        table.set_border_style(config.table_style);
        table << fort::header << "Point" << "Sign sequence" << "Sign variations" << fort::endr;
        for (size_t i = 0; i < sign_sequences.size(); i++) {
            table << (range[i] >= 0 ? " " + range[i].str() : range[i].str())
                  << SignSequenceToString(sign_sequences[i])
                  << sign_variations[i] << fort::endr;
        }
        std::cout << table.c_str();
    };

    const auto print_roots_table = [&]() {
        fort::char_table table;
        table.set_border_style(config.table_style);
        table << fort::header
              << "Number of real roots of p0(x) on interval (" + range[0].str() + ',' + range[1].str() + ']'
              << fort::endr;
        table << ("V(" + range[0].str() + ") - V(" + range[1].str() + ") = " + number_of_roots.str()) << fort::endr;
        std::cout << table.c_str();
    };

    if (vm.count("verbose")) {
        print_sturm_sequence();
        std::cout << std::endl;

        print_sign_table();
        std::cout << std::endl;

        print_roots_table();
    } else {
        std::cout << number_of_roots << std::endl;
    }
}

#endif