/**
 * @file   AtomicParameters.h
 *
 * @date   Jan 27, 2016
 * @author Thomas Dresselhaus
 *
 * @brief This file collects constant data for the elements (or isotopes) in the periodic table
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef ATOMICPARAMETERS_H_
#define ATOMICPARAMETERS_H_

/* Further includes */
/* Include Serenity Internal Headers */
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <map>
#include <vector>

namespace Serenity {

/************************************************************/
/* Unique data for elements, i.e. the same for all isotopes */
/************************************************************/

/* The periodic table in our program ends after Uuo, i.e. after the 7th period. */
constexpr unsigned int N_ELEMENTS_IN_PERIODIC_TABLE = 118;

/* The element symbols for the whole periodic table, the index matches the atomic number */
constexpr std::array<const char* const, N_ELEMENTS_IN_PERIODIC_TABLE + 1> PERIODIC_TABLE_OF_ELEMENTS{
    {// Dummy entry to make the array 1-based and make the index match the atomic number.
     // clang-format off
  "",
/* Group:
 * 1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 */
  "H" ,                                                                                "He",
  "Li","Be",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
  "Na","Mg",                                                  "Al","Si","P" ,"S" ,"Cl","Ar",
  "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
  "Cs","Ba","La",
                 "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",      // Lanthanoids
                 "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
  "Fr","Ra","Ac",
                 "Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",      // Actinoids
                 "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"}};
// clang-format on
/*
 * Each element has in its atomic form a certain, fixed configuration of the electron shell. This
 * is used e.g. in the AtomicScfGuessCalculator. Although the electron configuration looks very
 * systematic at first glance, there are quite a lot of exceptions caused by shielding effects of
 * inner electrons. Thus, we completely tabulate them to avoid any confusion.
 * These standard configurations for the elements have been taken from the wikipedia pages for each
 * particular element. (en.wikipedia.org, 26.01.2016)
 * Usage: For each element we list for each period how many electrons are in a shell with a certain
 * angular momentum.
 */
const std::vector<std::map<ANGULAR_QUANTUM_NUMBER, unsigned int>> ELEMENT_OCCUPATIONS[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and thus make the index match the atomic numbers
    {},
    // H
    {{{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // He
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Li
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Be
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // B
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // C
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // N
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // O
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // F
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
    // Ne
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}}},
    // Na
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Mg
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Al
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Si
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // P
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // S
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // Cl
    {{{ANGULAR_QUANTUM_NUMBER::s, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
// Ar
/* To shorten the entries below we put the argon configuration into a macro */
#define ARGON_OCCUPATION                                                \
  { ANGULAR_QUANTUM_NUMBER::s, 2 }                                      \
  }                                                                     \
  , {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}}, { \
    {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{ARGON_OCCUPATION}},
    // K
    {{ARGON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Ca
    {{ARGON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Sc
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 1}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ti
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // V
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 3}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Cr
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 5}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Mn
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 5}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Fe
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 6}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Co
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 7}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ni, an alternative occupation is [Ar] 3d9 4s1
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 8}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Cu
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Zn
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ga
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Ge
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // As
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // Se
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // Br
    {{ARGON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
// Kr
#define KRYPTON_OCCUPATION                            \
  ARGON_OCCUPATION, { ANGULAR_QUANTUM_NUMBER::d, 10 } \
  }                                                   \
  , {                                                 \
    {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{KRYPTON_OCCUPATION}},
    // Rb
    {{KRYPTON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Sr
    {{KRYPTON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Y
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 1}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Zr
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Nb
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 4}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Mo
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 5}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Tc
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 5}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ru
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 7}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Rh
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 8}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Pd
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}},
    // Ag
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Cd
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // In
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Zn
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // Sb
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // Te
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // I
    {{KRYPTON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 10}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
// Xe
#define XENON_OCCUPATION                                \
  KRYPTON_OCCUPATION, { ANGULAR_QUANTUM_NUMBER::d, 10 } \
  }                                                     \
  , {                                                   \
    {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{XENON_OCCUPATION}},
    // Cs
    {{XENON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Ba
    {{XENON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ln
    {{XENON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 1}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
// Ce
/*
 * In the lanthanoids the 4f shell is filled. To shorten the writing for that case we split up
 * the above macro into two parts
 */
#define XENON_OCC_PERIOD_1_TO_4 \
  KRYPTON_OCCUPATION, { ANGULAR_QUANTUM_NUMBER::d, 10 }
#define XENON_OCC_PERIOD_5 \
  {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 1}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Pr
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 3}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Nd
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 4}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Pm
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 5}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Sm
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 6}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Eu
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 7}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Gd
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 7}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Tb
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 9}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Dy
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 10}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ho
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 11}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Er
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 12}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Tm
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 13}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Yb
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}}, {XENON_OCC_PERIOD_5}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Lu
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Hf
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ta
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 3}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // W
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 4}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Re
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 5}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Os
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ir
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 7}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Pt
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 9}},
     {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Au
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Hg
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Tl
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Pb
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // Bi
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // Po
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // At
    {{XENON_OCC_PERIOD_1_TO_4, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
// Rn
#define RADON_OCCUPATION                                     \
  XENON_OCC_PERIOD_1_TO_4, { ANGULAR_QUANTUM_NUMBER::f, 14 } \
  }                                                          \
  , {XENON_OCC_PERIOD_5, {ANGULAR_QUANTUM_NUMBER::d, 10}}, { \
    {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{RADON_OCCUPATION}},
    // Fr
    {{RADON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 1}}},
    // Ra
    {{RADON_OCCUPATION}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ac
    {{RADON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 1}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Th
    {{RADON_OCCUPATION, {ANGULAR_QUANTUM_NUMBER::d, 2}}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
// Pa
#define RADON_OCC_PERIOD_1_TO_5 \
  XENON_OCCUPATION, { ANGULAR_QUANTUM_NUMBER::d, 10 }
#define RADON_OCC_PERIOD_6 \
  {ANGULAR_QUANTUM_NUMBER::s, 2}, { ANGULAR_QUANTUM_NUMBER::p, 6 }
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 2}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // U
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 3}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Np
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 4}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Pu
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 6}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Am
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 7}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Cm
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 7}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 1}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Bk
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 9}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Cf
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 10}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Es
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 11}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Fm
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 12}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Md
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 13}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // No
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}}, {RADON_OCC_PERIOD_6}, {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Lr
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Rf
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 2}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Db
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 3}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Sg
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 4}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Bh
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 5}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Hs
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 6}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Mt
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 7}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Ds
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 8}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Rg
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 9}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Cn
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}}},
    // Uut
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 1}}},
    // Fl
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 2}}},
    // Uup
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 3}}},
    // Lv
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 4}}},
    // Uus
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 5}}},
    // Uuo
    {{RADON_OCC_PERIOD_1_TO_5, {ANGULAR_QUANTUM_NUMBER::f, 14}},
     {RADON_OCC_PERIOD_6, {ANGULAR_QUANTUM_NUMBER::d, 10}},
     {{ANGULAR_QUANTUM_NUMBER::s, 2}, {ANGULAR_QUANTUM_NUMBER::p, 6}}},
};

#undef ARGON_OCCUPATION
#undef KRYPTON_OCCUPATION
#undef XENON_OCCUPATION
#undef XENON_OCC_PERIOD_1_TO_4
#undef XENON_OCC_PERIOD_5
#undef RADON_OCCUPATION
#undef RADON_OCC_PERIOD_1_TO_5
#undef RADON_OCC_PERIOD_6

/*
 * The data for these empirical sizes of the elements also listed on the wikipedia page
 * "Atomic radii of the elements (data page)" (en.wikipedia.org, 26.01.2016).
 * The values are given in an array but written in a style to represent the periodic table.
 * Data are in Angstrom.
 * In cases also that number is not available, a dummy entry of -1.0 is used which is properly
 * treated in the AtomType.h. The values are e.g. used when determining the grid accuracy.
 *
 * Original References:
 * J. Chem. Phys. 41, 3199 (main set)
 * J. Chem. Phys. 47, 1300 (use if no value in main set available)
 */
constexpr double ELEMENTAL_BRAGG_SLATER_RADII[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    // clang-format off
 -1.0,
/* Group:
 * 1 ,  2  ,  3  ,  4  ,  5  ,  6  ,  7  ,  8  ,  9  , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18 */
 0.25,                                                                                                 1.20,
 1.45, 1.05,                                                             0.85, 0.70, 0.65, 0.60, 0.50, 0.38,
 1.80, 1.50,                                                             1.25, 1.10, 1.00, 1.00, 1.00, 0.71,
 2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 0.88,
 2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.08,
 2.60, 2.15,
             1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, // Lanthanoids
                   1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.80, 1.60, 1.90, -1.0, 1.20,
 -1.0, 2.15,
             1.95, 1.80, 1.80, 1.75, 1.75, 1.75, 1.75, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, // Actinoids
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
};
// clang-format on
/*
 * The data for these empirical sizes of the elements are also listed on the wikipedia page
 * "Atomic radii of the elements (data page)" (en.wikipedia.org, 27.09.2016).
 * The values are given in an array but written in a style to represent the periodic table.
 * Data are in Angstrom.
 * In cases also that number is not available, a dummy entry of -1.0 is used which is properly
 * treated in the AtomType.h.
 *
 * Original References:
 *  J. Phys. Chem. 68, 441
 *  J. Phys. Chem. 2009, 113, 5806.
 */
constexpr double ELEMENTAL_VAN_DER_WAALS_RADII[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    // clang-format off
 -1.0,
/* Group:
 * 1 ,  2  ,  3  ,  4  ,  5  ,  6  ,  7  ,  8  ,  9  , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18 */
 1.20,                                                                                                 1.40,
 1.82, 1.53,                                                             1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
 2.27, 1.73,                                                             1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
 2.75, 2.31, 2.11, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.63, 1.40, 1.39, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02,
 3.03, 2.49, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.63, 1.72, 1.58, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16,
 3.43, 2.68,
             -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, // Lanthanoids
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.75, 1.66, 1.55, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20,
 3.48, 2.83,
             -1.0, -1.0, -1.0, 1.86, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, // Actinoids
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
};
// clang-format on
/*
 * UFF radii.\n
 * Reference:\n
 *   J. Am. Chem. Soc. 1992, 114, 25
 */
constexpr double ELEMENTAL_UFF_RADII[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    // clang-format off
 -1.0,
/* Group:
 *   1 ,    2  ,    3  ,    4   ,    5  ,    6  ,    7  ,    8  ,    9  ,   10  ,   11  ,   12  ,   13  ,   14  ,   15  ,   16  ,   17  ,   18 */
 1.4430,                                                                                                                                    1.81,
 1.2255, 1.3725,                                                                                  2.0515, 1.9255,   1.83, 1.75  , 1.682 , 1.6215,
 1.4915, 1.5105,                                                                                  2.2495, 2.1475, 2.0735, 2.0175, 1.9735,  1.934,
 1.9060, 1.6995, 1.6475, 1.5875, 1.5720, 1.5115 , 1.4805, 1.4560, 1.4360, 1.4170, 1.7475, 1.3815, 2.1915, 2.14  , 2.115 , 2.1025, 2.0945, 2.0705,
 2.0570, 1.8205, 1.6725, 1.5620, 1.5825, 1.526  , 1.499 , 1.4815, 1.4645, 1.4495, 1.5740, 1.4240, 2.2315, 2.1960, 2.2100, 2.2350, 2.25  , 2.2020,
 2.2585, 1.8515,
             1.7610, 1.7780, 1.8030, 1.7875, 1.7735, 1.7600, 1.7465, 1.6840, 1.7255, 1.7140, 1.7045, 1.6955, 1.6870, 1.6775, 1.8200, // Lanthanoids
                         1.5705, 1.850 , 1.5345, 1.4770, 1.5600 , 1.4200, 1.3770, 1.6465, 1.3525, 2.1735, 2.1485, 2.1850, 2.3545, 2.3750, 2.3825,
 2.4500, 1.8385,
             1.7390, 1.6980, 1.7120, 1.6975, 1.7120, 1.7120, 1.6905, 1.6630, 1.6695, 1.6565, 1.6495, 1.6430, 1.6370, 1.6240, 1.6180, // Actinoids
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
};
// clang-format on
/*
 * Number of core electrons for each element.\n
 * Adopted from the ORCA input documentation:\n
 *   https://sites.google.com/site/orcainputlibrary/frozen-core-calculations
 */
constexpr unsigned int NUMBER_OF_CORE_ELECTRONS[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    // clang-format off
    0,
/* Group:
 *   1 ,    2  ,    3  ,    4   ,    5  ,    6  ,    7  ,    8  ,    9  ,   10  ,   11  ,   12  ,   13  ,   14  ,   15  ,   16  ,   17  ,   18 */
     0,                                                                                                                                      0,
     0,   0,                                                                                         2,      2,     2,       2,      2,      2,
     2,   2,                                                                                        10,     10,    10,      10,     10,     10,
    10,  10,       10,     10,      10,     10,     10,     10,      10,     10,    10,      10,    18,     18,    18,      18,     18,     18,
    18,  18,       28,     28,      28,     28,     28,     28,      28,     28,    28,      28,    36,     36,    36,      36,     36,     36,
    36,  36,
    36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, // Lanthanoids
                   46,     46,      46,     46,     46,     46,      46,     46,    46,      46,    68,     68,     68,     68,     68,     68,
    68, 68,
    68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, // Actinoids
                   68,    100,     100,    100,    100,    100,     100,    100,    100,    100,    100,    100,    100,    100,    100,   100
};
// clang-format on

/**
 * Chemical hardnesses: d^2E/ dN^2 in au.
 * Source: International Journal of Quantum Chemistry, Vol 110, 1206 –1213 (2010) via sTDA GitHub page Oct 2021.
 * Note: These values are from the paper but multiplied by two since their definition is 1/2 d^2E/ dN^2, however
 * (IP-EA) = d^2 E/dN^2. See also on GitHub. stda/stda.f/line 1856.
 */
constexpr double CHEMICAL_HARDNESS[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    -1,          0.472592880, 0.922033910, 0.174528880, 0.257007330, 0.339490860, 0.421954120, 0.504381930, 0.586918630,
    0.669313510, 0.751916070, 0.179641050, 0.221572760, 0.263485780, 0.305396450, 0.347340140, 0.389247250, 0.431156700,
    0.473082690, 0.171054690, 0.202762440, 0.210073220, 0.217396470, 0.224710390, 0.232015010, 0.239339690, 0.246656380,
    0.253982550, 0.261288630, 0.268594760, 0.275925650, 0.307629990, 0.339315800, 0.372359850, 0.402735490, 0.434457760,
    0.466117080, 0.155850790, 0.186493240, 0.193562100, 0.200633110, 0.207705220, 0.214772540, 0.221846140, 0.228918720,
    0.235986210, 0.243056120, 0.250130180, 0.257199370, 0.287847800, 0.318486730, 0.349124310, 0.379765930, 0.410408080,
    0.441057770, 0.050193320, 0.067625700, 0.085044450, 0.102477360, 0.119911050, 0.137327720, 0.154762970, 0.172182650,
    0.189612880, 0.207047600, 0.224467520, 0.241896450, 0.259325030, 0.276760940, 0.294182310, 0.311595870, 0.329022740,
    0.345922980, 0.363880480, 0.381305860, 0.398774760, 0.416142980, 0.433645100, 0.451040140, 0.468489860, 0.485845500,
    0.125267300, 0.142686770, 0.160116150, 0.177558890, 0.194975570, 0.212407780, 0.072635250, 0.094221580, 0.099202950,
    0.104186210, 0.142356330, 0.163942940, 0.185519410, 0.223701390, -1.0,        -1.0,        -1.0,        -1.0,
    -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,
    -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,        -1.0,
    -1.0,        -1.0};
//clang-format on
/**
 * Covalent radii in Angstrom. Given radii are fitted for single-bond lengths.
 * Reference: Chem. Eur. J. 2008, 15, 186-197; https://doi.org/10.1002/chem.200800987
 */
constexpr double COVALENT_RADII[N_ELEMENTS_IN_PERIODIC_TABLE + 1] = {
    // Dummy entry to make the array 1-based and make the index match the atomic number.
    //clang-format off
    -1,
    /* Group:
     *   1 ,    2  ,    3  ,    4   ,    5  ,    6  ,    7  ,    8  ,    9  ,   10  ,   11  ,   12  ,   13  ,   14  , 15
     * ,   16  ,   17  ,   18 */
    0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67, 1.55, 1.39, 1.26, 1.16, 1.11, 1.03, 0.99, 0.96, 1.96,
    1.71, 1.48, 1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.10, 1.12, 1.18, 1.24, 1.21, 1.21, 1.16, 1.14, 1.17, 2.10, 1.85,
    1.63, 1.54, 1.47, 1.38, 1.28, 1.25, 1.25, 1.20, 1.28, 1.36, 1.42, 1.40, 1.40, 1.36, 1.33, 1.31, 2.32, 1.96, 1.62,
    1.52, 1.46, 1.37, 1.31, 1.29, 1.22, 1.23, 1.24, 1.33, 1.44, 1.44, 1.51, 1.45, 1.47, 1.42, 2.23, 2.01, 1.03, 1.04,
    1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.21, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18,
    // Lanthanoids
    1.80, 1.63, 1.76, 1.74, 1.73, 1.72, 1.68, 1.69, 1.68, 1.67, 1.66, 1.65, 1.64, 1.70,
    // Actinoids
    1.86, 1.75, 1.69, 1.70, 1.71, 1.72, 1.66, 1.66, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76};

} /* namespace Serenity */

#endif /* ATOMICPARAMETERS_H_ */
