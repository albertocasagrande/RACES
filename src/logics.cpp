/**
 * @file logics.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement a logic
 * @version 0.1
 * @date 2024-03-18
 *
 * @copyright Copyright (c) 2023-2024
 *
 * MIT License
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

#include "logics.hpp"

struct Context;

namespace Races
{

namespace Mutants
{

namespace Logics
{

Variable::Variable(const SpeciesId& species_id, const std::string& name):
    species_id(species_id), name(name)
{}

Expression::Expression(const Variable& variable):
    type(Type::VARIABLE), variable(variable), lhs(nullptr), rhs(nullptr)
{}

Expression::Expression(const size_t& value):
    type(Type::VALUE), value(value), lhs(nullptr), rhs(nullptr)
{}

Expression::Expression(const Type& type, const Expression& lhs, const Expression& rhs):
    type(type), lhs(std::make_shared<Expression>(lhs)),
    rhs(std::make_shared<Expression>(rhs))
{}

Expression::Expression(const Type& type, Expression&& lhs, Expression&& rhs):
    type(type), lhs(std::make_shared<Expression>(std::move(lhs))),
    rhs(std::make_shared<Expression>(std::move(rhs)))
{}

std::ostream& operator<<(std::ostream& os, const Expression& expression)
{
    switch(expression.type) {
        case Expression::Type::VALUE:
            os << expression.value;
            break;
        case Expression::Type::VARIABLE:
            os << expression.variable;
            break;
        case Expression::Type::SUM:
            os << *expression.lhs << "+" << *expression.rhs;
            break;
        case Expression::Type::SUB:
            os << *expression.lhs << "-" << *expression.rhs;
            break;
        case Expression::Type::MULT:
            os << *expression.lhs << "*" << *expression.rhs;
            break;
        default:
            throw std::runtime_error("Unknown expression type code " 
                                     + std::to_string(static_cast<unsigned int>(expression.type)));
    }

    return os;
}

Relation::Relation(const Type& type, Expression&& lhs, Expression&& rhs):
    type(type), lhs(lhs), rhs(rhs)
{}

Relation::Relation(const Type& type, const Expression& lhs, const Expression& rhs):
    type(type), lhs(lhs), rhs(rhs)
{}

std::ostream& operator<<(std::ostream& os, const Relation& relation)
{
    os << relation.lhs;
    switch(relation.type) {
        case Relation::Type::GT:
            os << ">";
            break;
        case Relation::Type::GE:
            os << ">=";
            break;
        case Relation::Type::EQ:
            os << "==";
            break;
        case Relation::Type::NE:
            os << "!=";
            break;
        default:
            throw std::runtime_error("Unknown relation type code " 
                                     + std::to_string(static_cast<unsigned int>(relation.type)));
    }
    os << relation.rhs;

    return os;
}


Formula::Formula():
    rel(nullptr), lhs(nullptr), rhs(nullptr)
{}

Formula::Formula(const Type& type, const Formula& lhs, const Formula& rhs):
    type(type), rel(nullptr), lhs(std::make_shared<Formula>(lhs)),
    rhs(std::make_shared<Formula>(rhs))
{
    if (type == Type::RELATION) {
        throw std::domain_error("Building a formula with type relation "
                                "requires a relation as the parameter.");
    }
    if (type == Type::NEG) {
        throw std::domain_error("Building a formula with type negation "
                                "requires a single formula as the parameter.");
    }
}

Formula::Formula(const Type& type, const Formula& subformula):
    type(type), rel(nullptr), lhs(std::make_shared<Formula>(subformula)),
    rhs(nullptr)
{
    if (type == Type::RELATION) {
        throw std::domain_error("Building a formula with type relation "
                                "requires a relation as the parameter.");
    }
    if (type != Type::NEG) {
        throw std::domain_error("Building a formula with a type different "
                                "from negation requires two formulas.");
    }
}

Formula::Formula(const Relation& relation):
    type(Type::RELATION), rel(std::make_shared<Relation>(relation)),
    lhs(nullptr), rhs(nullptr)
{}

Formula::Formula(const Type& type, Formula&& lhs, Formula&& rhs):
    type(type), rel(nullptr), lhs(std::make_shared<Formula>(std::move(lhs))),
    rhs(std::make_shared<Formula>(std::move(rhs)))
{
    if (type == Type::RELATION) {
        throw std::domain_error("Building a formula with type relation "
                                "requires a relation as the parameter.");
    }
    if (type == Type::NEG) {
        throw std::domain_error("Building a formula with type negation "
                                "requires a single formula as the parameter.");
    }
}

Formula::Formula(const Type& type, Formula&& subformula):
    type(Type::NEG), rel(nullptr), lhs(std::make_shared<Formula>(std::move(subformula))),
    rhs(nullptr)
{
    if (type == Type::RELATION) {
        throw std::domain_error("Building a formula with type relation "
                                "requires a relation as the parameter.");
    }
    if (type != Type::NEG) {
        throw std::domain_error("Building a formula with a type different "
                                "from negation requires two formulas.");
    }
}

Formula::Formula(Relation&& relation):
    type(Type::RELATION), rel(std::make_shared<Relation>(std::move(relation))),
    lhs(nullptr), rhs(nullptr)
{}

std::ostream& operator<<(std::ostream& os, const Formula& formula)
{
    switch(formula.type) {
        case Formula::Type::RELATION:
            os << *(formula.rel);
            break;
        case Formula::Type::NEG:
            os << "!(" << *(formula.lhs) << ")";
            break;
        case Formula::Type::AND:
            os << "(" << *(formula.lhs) <<  "&&" << *(formula.rhs) << ")";
            break;
        case Formula::Type::OR:
            os << "(" << *(formula.lhs) <<  "||" << *(formula.rhs) << ")";
            break;
        default:
            throw std::runtime_error("Unknown formula type code " 
                                        + std::to_string(static_cast<unsigned int>(formula.type)));
    }

    return os;
}


}   // Logics

}   // Mutants

}   // Races