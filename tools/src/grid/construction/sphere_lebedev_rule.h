/**
 * @file   sphere_lebedev_rule.h
 *
 * @date  Jun 26, 2015 (last cleaning)
 * @author Thomas Dresselhaus, last cleaning Jan Unsleber
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
 *
 *  This code was taken and modified from:
 *  http://people.sc.fsu.edu/~jburkardt/f_src/sphere_lebedev_rule/sphere_lebedev_rule.html
 *  and
 *  http://people.sc.fsu.edu/~jburkardt/f_src/sphere_lebedev_rule/sphere_lebedev_rule.f90
 *  (last access date: Aug 24, 2017)
 *
 *  The original Code: "SPHERE_LEBEDEV_RULE" was/is licensed under the:
 *  GNU Lesser General Public License (L-GPL)
 *
 */

#ifndef SPHERE_LEBEDEV_RULE_H_
#define SPHERE_LEBEDEV_RULE_H_

namespace Serenity {

int gen_oh(int code, double a, double b, double v, double* x, double* y, double* z, double* w);
void ld0006(double* x, double* y, double* z, double* w);
void ld0014(double* x, double* y, double* z, double* w);
void ld0026(double* x, double* y, double* z, double* w);
void ld0038(double* x, double* y, double* z, double* w);
void ld0050(double* x, double* y, double* z, double* w);
void ld0074(double* x, double* y, double* z, double* w);
void ld0086(double* x, double* y, double* z, double* w);
void ld0110(double* x, double* y, double* z, double* w);
void ld0146(double* x, double* y, double* z, double* w);
void ld0170(double* x, double* y, double* z, double* w);
void ld0194(double* x, double* y, double* z, double* w);
void ld0230(double* x, double* y, double* z, double* w);
void ld0266(double* x, double* y, double* z, double* w);
void ld0302(double* x, double* y, double* z, double* w);
void ld0350(double* x, double* y, double* z, double* w);
void ld0434(double* x, double* y, double* z, double* w);
void ld0590(double* x, double* y, double* z, double* w);
void ld0770(double* x, double* y, double* z, double* w);
void ld0974(double* x, double* y, double* z, double* w);
void ld1202(double* x, double* y, double* z, double* w);
void ld1454(double* x, double* y, double* z, double* w);
void ld1730(double* x, double* y, double* z, double* w);
void ld2030(double* x, double* y, double* z, double* w);
void ld2354(double* x, double* y, double* z, double* w);
void ld2702(double* x, double* y, double* z, double* w);
void ld3074(double* x, double* y, double* z, double* w);
void ld3470(double* x, double* y, double* z, double* w);
void ld3890(double* x, double* y, double* z, double* w);
void ld4334(double* x, double* y, double* z, double* w);
void ld4802(double* x, double* y, double* z, double* w);
void ld5294(double* x, double* y, double* z, double* w);
void ld5810(double* x, double* y, double* z, double* w);

} /*namespace Serenity*/
#endif /* SPHERE_LEBEDEV_RULE_H_ */
