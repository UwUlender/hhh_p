#pragma once
#include <string>
#include <vector>
#include <map>
#include <stack>
#include <stdexcept>
#include <cmath>

namespace HydroPlas {

enum class TokenType { Number, Variable, Operator, Function, LeftParen, RightParen };

struct Token {
    TokenType type;
    std::string value;
    double num_val = 0.0;
};

class EquationEvaluator {
public:
    EquationEvaluator() = default;
    
    // Parse the equation string (infix) into RPN
    void parse(const std::string& equation);
    
    // Evaluate the parsed RPN using the provided variable map
    double evaluate(const std::map<std::string, double>& variables) const;
    
    // Set constants that don't change
    void set_constant(const std::string& name, double value);

private:
    std::vector<Token> rpn_queue_;
    std::map<std::string, double> constants_;
    
    int get_precedence(const std::string& op) const;
    bool is_right_assoc(const std::string& op) const;
    bool is_operator(char c) const;
};

} // namespace HydroPlas
