#include "EquationEvaluator.hpp"
#include <iostream>
#include <sstream>
#include <cctype>
#include <algorithm>

namespace HydroPlas {

void EquationEvaluator::set_constant(const std::string& name, double value) {
    constants_[name] = value;
}

int EquationEvaluator::get_precedence(const std::string& op) const {
    if (op == "+" || op == "-") return 1;
    if (op == "*" || op == "/") return 2;
    if (op == "^") return 3;
    return 0;
}

bool EquationEvaluator::is_right_assoc(const std::string& op) const {
    return op == "^";
}

bool EquationEvaluator::is_operator(char c) const {
    return c == '+' || c == '-' || c == '*' || c == '/' || c == '^';
}

void EquationEvaluator::parse(const std::string& equation) {
    rpn_queue_.clear();
    std::stack<Token> op_stack;
    
    // Remove outer braces {} if present (YAML specific)
    std::string eq = equation;
    if (eq.front() == '{' && eq.back() == '}') {
        eq = eq.substr(1, eq.length() - 2);
    }
    
    size_t i = 0;
    while (i < eq.length()) {
        char c = eq[i];
        
        if (std::isspace(c)) {
            i++;
            continue;
        }
        
        // Number
        if (std::isdigit(c) || c == '.') {
            std::string num_str;
            while (i < eq.length() && (std::isdigit(eq[i]) || eq[i] == '.' || eq[i] == 'e' || eq[i] == 'E')) {
                // Handle scientific notation + or - (e.g., 1.0e-3)
                if ((eq[i] == 'e' || eq[i] == 'E') && i + 1 < eq.length() && (eq[i+1] == '+' || eq[i+1] == '-')) {
                    num_str += eq[i];
                    i++;
                    num_str += eq[i];
                    i++;
                    continue;
                }
                num_str += eq[i];
                i++;
            }
            try {
                Token t;
                t.type = TokenType::Number;
                t.num_val = std::stod(num_str);
                rpn_queue_.push_back(t);
            } catch (...) {
                throw std::runtime_error("Invalid number format: " + num_str);
            }
        }
        // Variable or Function
        else if (std::isalpha(c) || c == '_') {
            std::string name;
            while (i < eq.length() && (std::isalnum(eq[i]) || eq[i] == '_')) {
                name += eq[i];
                i++;
            }
            
            // Check if function
            if (name == "exp" || name == "log" || name == "sin" || name == "cos" || name == "sqrt") {
                Token t;
                t.type = TokenType::Function;
                t.value = name;
                op_stack.push(t);
            } else {
                Token t;
                t.type = TokenType::Variable;
                t.value = name;
                rpn_queue_.push_back(t);
            }
        }
        // Operator
        else if (is_operator(c)) {
            // Check for unary minus
            // It's unary if: start of string, or previous token was Operator or LeftParen
            bool is_unary = false;
            if (c == '-') {
                if (rpn_queue_.empty()) is_unary = true; // Start of expression
                // Need to check what was the last thing added to output? No, depends on infix structure.
                // We don't track infix history easily here.
                // Simplified check: If we just read an operator or paren, it's unary.
                // This parser is simple. Let's assume standard infix.
                // To support unary minus properly in Shunting Yard is tricky.
                // Hack: Replace "-x" with "(0-x)" or treat as special operator "u-".
                // Given the example "exp((-1.0e-3)...)", the minus is part of the number parsing above?
                // "(-1.0e-3)" -> If '(' is processed, next is '-'. 
                // If I see '-' and previous was not a number/variable/RightParen, it's unary.
            }

            // For now, assume space-separated or standard.
            // Let's implement unary minus as a special operator '#'
            // But checking context is hard without previous token type.
            // Let's assume the user puts 0-x or -1*x if needed, OR simple unary support.
            // Example: "(-1.0e-3)" is parsed as number "-1.0e-3" if the loop catches it?
            // My number parser only starts with digit or dot. It doesn't start with '-'.
            // So '-' will be seen here.
            
            // Correct logic for unary minus:
            // If we are at start (i==0/skipped spaces) or prev char was '(' or operator.
            // We need to look back at 'last_token_type'.
            // Implementation detail: I'll skip complex unary minus for now and rely on "0 - ..." or "-1 * ..." 
            // BUT the example has "exp((-1.0e-3)...)". 
            // My number parser should handle signed numbers if they start with -.
            // No, standard `stod` handles it, but my tokenizer splits on '-'.
            // Workaround: If c is '-' and next char is digit, treat as number start?
            if (c == '-' && i+1 < eq.length() && (std::isdigit(eq[i+1]) || eq[i+1] == '.')) {
                 // It is a negative number!
                 std::string num_str = "-";
                 i++; // Consume '-'
                 while (i < eq.length() && (std::isdigit(eq[i]) || eq[i] == '.' || eq[i] == 'e' || eq[i] == 'E')) {
                    if ((eq[i] == 'e' || eq[i] == 'E') && i + 1 < eq.length() && (eq[i+1] == '+' || eq[i+1] == '-')) {
                        num_str += eq[i]; i++; num_str += eq[i]; i++; continue;
                    }
                    num_str += eq[i]; i++;
                 }
                 Token t; t.type = TokenType::Number; t.num_val = std::stod(num_str);
                 rpn_queue_.push_back(t);
                 continue;
            }

            std::string op(1, c);
            Token t; t.type = TokenType::Operator; t.value = op;
            
            while (!op_stack.empty() && op_stack.top().type != TokenType::LeftParen) {
                const Token& top = op_stack.top();
                if (top.type == TokenType::Function ||
                    get_precedence(top.value) > get_precedence(op) ||
                    (get_precedence(top.value) == get_precedence(op) && !is_right_assoc(op))) {
                    rpn_queue_.push_back(top);
                    op_stack.pop();
                } else {
                    break;
                }
            }
            op_stack.push(t);
            i++;
        }
        // Parentheses
        else if (c == '(') {
            Token t; t.type = TokenType::LeftParen;
            op_stack.push(t);
            i++;
        }
        else if (c == ')') {
            while (!op_stack.empty() && op_stack.top().type != TokenType::LeftParen) {
                rpn_queue_.push_back(op_stack.top());
                op_stack.pop();
            }
            if (!op_stack.empty()) op_stack.pop(); // Pop '('
            // Check if function
            if (!op_stack.empty() && op_stack.top().type == TokenType::Function) {
                rpn_queue_.push_back(op_stack.top());
                op_stack.pop();
            }
            i++;
        }
        else {
            i++; // Skip unknown chars
        }
    }
    
    while (!op_stack.empty()) {
        rpn_queue_.push_back(op_stack.top());
        op_stack.pop();
    }
}

double EquationEvaluator::evaluate(const std::map<std::string, double>& variables) const {
    std::stack<double> val_stack;
    
    for (const auto& token : rpn_queue_) {
        if (token.type == TokenType::Number) {
            val_stack.push(token.num_val);
        }
        else if (token.type == TokenType::Variable) {
            // Check provided variables first
            auto it = variables.find(token.value);
            if (it != variables.end()) {
                val_stack.push(it->second);
            } else {
                // Check constants
                auto it_c = constants_.find(token.value);
                if (it_c != constants_.end()) {
                    val_stack.push(it_c->second);
                } else {
                    // Default to 0.0 or throw?
                    // throw std::runtime_error("Unknown variable: " + token.value);
                    val_stack.push(0.0);
                }
            }
        }
        else if (token.type == TokenType::Operator) {
            if (val_stack.size() < 2) throw std::runtime_error("Stack underflow for op " + token.value);
            double b = val_stack.top(); val_stack.pop();
            double a = val_stack.top(); val_stack.pop();
            
            if (token.value == "+") val_stack.push(a + b);
            else if (token.value == "-") val_stack.push(a - b);
            else if (token.value == "*") val_stack.push(a * b);
            else if (token.value == "/") val_stack.push(a / b);
            else if (token.value == "^") val_stack.push(std::pow(a, b));
        }
        else if (token.type == TokenType::Function) {
            if (val_stack.empty()) throw std::runtime_error("Stack underflow for func " + token.value);
            double a = val_stack.top(); val_stack.pop();
            
            if (token.value == "exp") val_stack.push(std::exp(a));
            else if (token.value == "log") val_stack.push(std::log(a));
            else if (token.value == "sin") val_stack.push(std::sin(a));
            else if (token.value == "cos") val_stack.push(std::cos(a));
            else if (token.value == "sqrt") val_stack.push(std::sqrt(a));
        }
    }
    
    return val_stack.empty() ? 0.0 : val_stack.top();
}

} // namespace HydroPlas
