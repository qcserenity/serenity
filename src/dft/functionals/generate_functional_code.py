#!/usr/bin/python3
# @file   generate_functional_code.py
#
# @date   Sep 3, 2020
# @author Jan Unsleber
# @copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n
#
# =============================================================
#
#  This file generates most of the code that is used to build
#  and resolve functionals in Serenity.
#
#  It requires the following files:
#   1. basic_functionals.dat
#   2. composite_functionals.dat
#   3. functional_definitions.dat
#
#  It generates the following files:
#   1. BasicFunctionals.h
#   2. BasicFunctionals.cpp
#   3. BasicFunctionals_tests.cpp
#   4. BasicFunctionals_python.pp
#   5. CompositeFunctionals.h
#   6. CompositeFunctionals.cpp
#   7. CompositeFunctionals_tests.cpp
#   8. CompositeFunctionals_python.pp
#
# =============================================================

functionals = []
with open('basic_functionals.dat', 'r') as f:
    for l in f:
        if l[0] != '#':
            functionals.append(l)
nFuncs = len(functionals)

with open('BasicFunctionals.h', 'w') as bf:
    bf.write('/**\n')
    bf.write(' * @file   BasicFunctionals.h\n')
    bf.write(' *\n')
    bf.write(' * @date   Sep 3, 2020\n')
    bf.write(' * @author Jan P. Unsleber\n')
    bf.write(' *\n')
    bf.write(' * IMPORTANT:\\n\n')
    bf.write(' * This file was automatically generated, please do not alter it.\n')
    bf.write(' * Any required changes should be made to the generating Python script\n')
    bf.write(' * which should be located close by.\n')
    bf.write(' *\n')
    bf.write(' * @copyright \\n\n')
    bf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    bf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    bf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    bf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    bf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    bf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    bf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    bf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    bf.write(' *  GNU General Public License for more details.\\n\\n\n')
    bf.write(' *  You should have received a copy of the GNU Lesser General\n')
    bf.write(' *  Public License along with Serenity.\n')
    bf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    bf.write(' */\n')
    bf.write('#ifndef BASICFUNCTIONALS_H_\n')
    bf.write('#define BASICFUNCTIONALS_H_\n\n')
    bf.write('/* Include Std and External Headers */\n')
    bf.write('#include <array>\n\n')
    bf.write('namespace Serenity {\n')
    bf.write('namespace BasicFunctionals {\n\n')
    bf.write('enum class PURPOSES {KINETIC, EXCHANGE, CORRELATION, EXCHANGE_CORRELATION, NONE};\n')
    bf.write('enum class CLASSES {NONE, LDA, GGA, META_GGA, MODELL};\n')
    bf.write('enum class IMPLEMENTATIONS {XCFUN, LIBXC, BOTH};\n\n')

    # Define BASIC_FUNCTIONALS
    bf.write('enum class BASIC_FUNCTIONALS {\n')
    for i, f in enumerate(functionals):
        data = f.split()
        bf.write('  {} = {:d}'.format(data[3], i))
        if i < nFuncs:
            bf.write(',\n')
        else:
            bf.write('\n')
    bf.write('};\n')

    bf.write('#ifdef SERENITY_USE_LIBXC\n')
    bf.write('int getLibXCAlias(BASIC_FUNCTIONALS& functional);\n')
    bf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    bf.write('#ifdef SERENITY_USE_XCFUN\n')
    bf.write('char* getXCFunAlias(BASIC_FUNCTIONALS& functional);\n')
    bf.write('#endif /* SERENITY_USE_XCFUN */\n\n')

    # Define Purposes
    def resolve_purpose(s):
        if s == 'K':
            return 'PURPOSES::KINETIC'
        if s == 'X':
            return 'PURPOSES::EXCHANGE'
        if s == 'C':
            return 'PURPOSES::CORRELATION'
        if s == 'XC':
            return 'PURPOSES::EXCHANGE_CORRELATION'
        if s == 'N':
            return 'PURPOSES::NONE'
        raise Exception("Wrong Purpose")

    bf.write('constexpr std::array<PURPOSES, '+str(nFuncs)+'> getPurpose{\n')
    for i, f in enumerate(functionals):
        data = f.split()
        bf.write('  {}  /* {:d} */'.format(resolve_purpose(data[1]), i))
        if i < nFuncs:
            bf.write(',\n')
        else:
            bf.write('\n')
    bf.write('};\n\n')

    # Define Class
    def resolve_class(s):
        if s == 'MOD':
            return 'CLASSES::MODELL'
        if s == 'MGGA':
            return 'CLASSES::META_GGA'
        return 'CLASSES::'+s

    bf.write('constexpr std::array<CLASSES,'+str(nFuncs)+'> getClass{\n')
    for i, f in enumerate(functionals):
        data = f.split()
        bf.write('  {} /* {:d} */'.format(resolve_class(data[0]), i))
        if i < nFuncs:
            bf.write(',\n')
        else:
            bf.write('\n')
    bf.write('};\n\n')

    # Define IMPLEMENTATIONS
    bf.write('constexpr std::array<IMPLEMENTATIONS, '+str(nFuncs)+'> getImplementation{\n')
    for i, f in enumerate(functionals):
        data = f.split()
        bf.write('  IMPLEMENTATIONS::{} /* {:d} */'.format(data[2], i))
        if i < nFuncs:
            bf.write(',\n')
        else:
            bf.write('\n')
    bf.write('};\n\n')
    bf.write('} /* namespace BasicFunctionals */\n')
    bf.write('} /* namespace Serenity */\n\n')
    bf.write('#endif /* BASICFUNCTIONALS_H_ */\n')

with open('BasicFunctionals.cpp', 'w') as bf:
    bf.write('/**\n')
    bf.write(' * @file   BasicFunctionals.cpp\n')
    bf.write(' *\n')
    bf.write(' * @date   Sep 3, 2020\n')
    bf.write(' * @author Jan P. Unsleber\n')
    bf.write(' *\n')
    bf.write(' * IMPORTANT:\\n\n')
    bf.write(' * This file was automatically generated, please do not alter it.\n')
    bf.write(' * Any required changes should be made to the generating Python script\n')
    bf.write(' * which should be located close by.\n')
    bf.write(' *\n')
    bf.write(' * @copyright \\n\n')
    bf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    bf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    bf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    bf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    bf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    bf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    bf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    bf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    bf.write(' *  GNU General Public License for more details.\\n\\n\n')
    bf.write(' *  You should have received a copy of the GNU Lesser General\n')
    bf.write(' *  Public License along with Serenity.\n')
    bf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    bf.write(' */\n')
    bf.write('/* Include Class Header*/\n')
    bf.write('#include "dft/functionals/BasicFunctionals.h"\n')
    bf.write('/* Include Serenity Internal Headers */\n')
    bf.write('#include "misc/SerenityError.h"\n')
    bf.write('#include "settings/DFTOptions.h"\n')
    bf.write('/* Include Std and External Headers */\n')
    bf.write('#include <string>\n')
    bf.write('#ifdef SERENITY_USE_LIBXC\n')
    bf.write('#include <xc_funcs.h>\n')
    bf.write('#endif /* SERENITY_USE_LIBXC */\n')
    bf.write('\n')
    bf.write('namespace Serenity {\n')
    bf.write('namespace BasicFunctionals {\n\n')
    bf.write('#ifdef SERENITY_USE_LIBXC\n')
    bf.write('int getLibXCAlias(BASIC_FUNCTIONALS& functional) {\n')
    bf.write('  int alias = -1;\n')
    bf.write('  switch (functional) {\n')
    for f in functionals:
        data = f.split()
        if data[4] == '-':
            continue
        bf.write('    case BASIC_FUNCTIONALS::{}:\n'.format(data[3]))
        bf.write('      alias = {};\n'.format(data[4]))
        bf.write('      break;\n')
    bf.write('    default:\n')
    bf.write('      std::string s;\n')
    bf.write('      Options::resolve<BasicFunctionals::BASIC_FUNCTIONALS>(s, functional);\n')
    bf.write('      throw SerenityError("Basic functional " + s + " unknown to LibXC.");\n')
    bf.write('      break;\n')
    bf.write('  }\n')
    bf.write('  return alias;\n')
    bf.write('}\n')
    bf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    bf.write('#ifdef SERENITY_USE_XCFUN\n')
    bf.write('char* getXCFunAlias(BASIC_FUNCTIONALS& functional) {\n')
    bf.write('  char* alias = (char*)"null";\n')
    bf.write('  switch (functional) {\n')
    for f in functionals:
        data = f.split()
        if data[5] == '-':
            continue
        bf.write('    case BASIC_FUNCTIONALS::{}:\n'.format(data[3]))
        bf.write('      alias = (char*)"{}";\n'.format(data[5]))
        bf.write('      break;\n')
    bf.write('    default:\n')
    bf.write('      std::string s;\n')
    bf.write('      Options::resolve<BasicFunctionals::BASIC_FUNCTIONALS>(s, functional);\n')
    bf.write('      throw SerenityError("Basic functional " + s + " unknown to XCFun.");\n')
    bf.write('      break;\n')
    bf.write('  }\n')
    bf.write('  return alias;\n')
    bf.write('}\n')
    bf.write('#endif /* SERENITY_USE_XCFUN */\n\n')
    bf.write('} /* namespace BasicFunctionals */\n')
    bf.write('namespace Options {\n')
    bf.write('template<>')
    bf.write('void resolve<BasicFunctionals::BASIC_FUNCTIONALS>(std::string& value, BasicFunctionals::BASIC_FUNCTIONALS& field) {')
    bf.write('static const std::map<std::string, BasicFunctionals::BASIC_FUNCTIONALS> m = {\n')
    first = True
    for f in functionals:
        data = f.split()
        if not first:
            bf.write(',\n')
        bf.write('    {{"{}", BasicFunctionals::BASIC_FUNCTIONALS::{}}}'.format(data[3], data[3]))
        first = False
    bf.write('};\n')
    bf.write('check(m, value, field);\n')
    bf.write('}\n')
    bf.write('} /* namespace Options */\n')
    bf.write('} /* namespace Serenity */\n\n')

functionals = []
with open('composite_functionals.dat', 'r') as f:
    for l in f:
        if l[0] != '#':
            functionals.append(l)
nFuncs = len(functionals)
basics = []
with open('basic_functionals.dat', 'r') as f:
    for l in f:
        if l[0] != '#':
            basics.append(l)
nBFuncs = len(basics)
definitions = []
with open('functional_definitions.dat', 'r') as f:
    for l in f:
        if l[0] != '#':
            definitions.append(l)

with open('CompositeFunctionals.h', 'w') as cf:
    cf.write('/**\n')
    cf.write(' * @file   CompositeFunctionals.h\n')
    cf.write(' *\n')
    cf.write(' * @date   Sep 3, 2020\n')
    cf.write(' * @author Jan P. Unsleber\n')
    cf.write(' *\n')
    cf.write(' * IMPORTANT:\\n\n')
    cf.write(' * This file was automatically generated, please do not alter it.\n')
    cf.write(' * Any required changes should be made to the generating Python script\n')
    cf.write(' * which should be located close by.\n')
    cf.write(' *\n')
    cf.write(' * @copyright \\n\n')
    cf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    cf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    cf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    cf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    cf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    cf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    cf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    cf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    cf.write(' *  GNU General Public License for more details.\\n\\n\n')
    cf.write(' *  You should have received a copy of the GNU Lesser General\n')
    cf.write(' *  Public License along with Serenity.\n')
    cf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    cf.write(' */\n')
    cf.write('#ifndef COMPOSITEFUNCTIONALS_H_\n')
    cf.write('#define COMPOSITEFUNCTIONALS_H_\n')
    cf.write('/* Include Std and External Headers */\n')
    cf.write('#include <array>\n\n')
    cf.write('\n')
    cf.write('namespace Serenity {\n\n')
    cf.write('class Functional;\n\n')
    cf.write('namespace CompositeFunctionals {\n\n')
    cf.write('enum class PURPOSES {KINETIC, EXCHANGE_CORRELATION, NONE};\n')
    cf.write('enum class CLASSES {NONE, LDA, GGA, META_GGA, MODELL};\n')
    cf.write('enum class IMPLEMENTATIONS {XCFUN, LIBXC, EITHER_OR, BOTH};\n\n')
    cf.write('enum class FUNCTIONALS {\n')
    for i, f in enumerate(functionals):
        data = f.split()
        cf.write('  {} = {:d}'.format(data[0], i))
        if i < nFuncs:
            cf.write(',\n')
        else:
            cf.write('\n')
    cf.write('};\n\n')
    cf.write('enum class XCFUNCTIONALS {\n')
    for i, f in enumerate(functionals):
        data = f.split()
        if not data[1] == 'XC' and not data[1] == 'NONE':
            continue
        cf.write('  {} = {:d}'.format(data[0], i))
        if i < nFuncs:
            cf.write(',\n')
        else:
            cf.write('\n')
    cf.write('};\n\n')
    cf.write('enum class KINFUNCTIONALS {\n')
    for i, f in enumerate(functionals):
        data = f.split()
        if not data[1] == 'K' and not data[1] == 'NONE':
            continue
        cf.write('  {} = {:d}'.format(data[0], i))
        if i < nFuncs:
            cf.write(',\n')
        else:
            cf.write('\n')
    cf.write('};\n\n')
    cf.write('Functional resolveFunctional(FUNCTIONALS functional);\n\n')
    cf.write('Functional resolveFunctional(XCFUNCTIONALS functional);\n\n')
    cf.write('Functional resolveFunctional(KINFUNCTIONALS functional);\n\n')
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('Functional resolveLibXC(FUNCTIONALS functional);\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('Functional resolveXCFun(FUNCTIONALS functional);\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n\n')
    # Define Purposes
    cf.write('constexpr std::array<PURPOSES, '+str(nFuncs)+'> getPurpose{\n')
    for i, f in enumerate(functionals):
        data = f.split()
        if data[1] == 'K':
            cf.write('  {}  /* {:d} */'.format('PURPOSES::KINETIC', i))
        elif data[1] == 'XC':
            cf.write('  {}  /* {:d} */'.format('PURPOSES::EXCHANGE_CORRELATION', i))
        elif data[1] == 'NONE':
            cf.write('  {}  /* {:d} */'.format('PURPOSES::NONE', i))
        else:
            raise RuntimeError('False purpose in Composites.')
        if i < nFuncs:
            cf.write(',\n')
        else:
            cf.write('\n')
    cf.write('};\n\n')
    # Define IMPLEMENTATIONS
    cf.write('constexpr std::array<IMPLEMENTATIONS, '+str(nFuncs)+'> getImplementation{\n')
    impl_tot = {}
    for d in definitions:
        split = d.split('|')
        f = split[0].split()[0]
        impl = split[0].split()[1]
        if f in impl_tot:
            impl_tot[f] = 'EITHER_OR'
        else:
            impl_tot[f] = impl
    for i, f in enumerate(functionals):
        data = f.split()
        cf.write('  IMPLEMENTATIONS::{} /* {:d} */'.format(impl_tot[data[0]], i))
        if i < nFuncs:
            cf.write(',\n')
        else:
            cf.write('\n')
    cf.write('};\n\n')
    cf.write('} /* namespace CompositeFunctionals */\n')
    cf.write('} /* namespace Serenity */\n\n')
    cf.write('#endif /* COMPOSITEFUNCTIONALS_H_ */\n')


with open('CompositeFunctionals.cpp', 'w') as cf:
    cf.write('/**\n')
    cf.write(' * @file   CompositeFunctionals.cpp\n')
    cf.write(' *\n')
    cf.write(' * @date   Sep 3, 2020\n')
    cf.write(' * @author Jan P. Unsleber\n')
    cf.write(' *\n')
    cf.write(' * IMPORTANT:\\n\n')
    cf.write(' * This file was automatically generated, please do not alter it.\n')
    cf.write(' * Any required changes should be made to the generating Python script\n')
    cf.write(' * which should be located close by.\n')
    cf.write(' *\n')
    cf.write(' * @copyright \\n\n')
    cf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    cf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    cf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    cf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    cf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    cf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    cf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    cf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    cf.write(' *  GNU General Public License for more details.\\n\\n\n')
    cf.write(' *  You should have received a copy of the GNU Lesser General\n')
    cf.write(' *  Public License along with Serenity.\n')
    cf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    cf.write(' */\n')
    cf.write('/* Include Class Header*/\n')
    cf.write('#include "dft/functionals/CompositeFunctionals.h"\n')
    cf.write('/* Include Serenity Internal Headers */\n')
    cf.write('#include "dft/Functional.h"\n')
    cf.write('#include "dft/functionals/BasicFunctionals.h"\n')
    cf.write('#include "misc/SerenityError.h"\n')
    cf.write('#include "settings/Options.h"\n')
    cf.write('/* Include Std and External Headers */\n')
    cf.write('#include <map>\n')
    cf.write('#include <string>\n')
    cf.write('#include <utility>\n')
    cf.write('#include <vector>\n')
    cf.write('\n')
    cf.write('namespace Serenity {\n')
    cf.write('namespace CompositeFunctionals {\n\n')
    cf.write('Functional resolveFunctional(FUNCTIONALS functional) {\n')
    cf.write('#if defined SERENITY_PREFER_XCFUN && defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC\n')
    cf.write('  try {\n')
    cf.write('    return resolveXCFun(functional);\n')
    cf.write('  }\n')
    cf.write('  catch (...) {\n')
    cf.write('    try {\n')
    cf.write('      return resolveLibXC(functional);\n')
    cf.write('    }\n')
    cf.write('    catch (...) {\n')
    cf.write('      throw SerenityError("You have requested a functional that neither XCFun nor LibXC can provide.");\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('#elif defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC\n')
    cf.write('  try {\n')
    cf.write('    return resolveLibXC(functional);\n')
    cf.write('  }\n')
    cf.write('  catch (...) {\n')
    cf.write('    try {\n')
    cf.write('      return resolveXCFun(functional);\n')
    cf.write('    }\n')
    cf.write('    catch (...) {\n')
    cf.write('      throw SerenityError("You have requested a functional that neither XCFun nor LibXC can provide.");\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('#elif defined SERENITY_USE_XCFUN\n')
    cf.write('  return resolveXCFun(functional);\n')
    cf.write('#else\n')
    cf.write('  return resolveLibXC(functional);\n')
    cf.write('#endif\n')
    cf.write('}\n\n')
    cf.write('Functional resolveFunctional(XCFUNCTIONALS functional) {\n')
    cf.write('  return resolveFunctional(FUNCTIONALS((int)functional));\n')
    cf.write('}\n\n')
    cf.write('Functional resolveFunctional(KINFUNCTIONALS functional) {\n')
    cf.write('  return resolveFunctional(FUNCTIONALS((int)functional));\n')
    cf.write('}\n\n')
    # Parse Functional definitions and check against other lists
    for d in definitions:
        split = d.split('|')
        f = split[0].split()[0]
        impl = split[0].split()[1]
        bfs = split[1].split()
        bfmix = split[2].split()
        factors = split[3].split()
        if len(bfs) != len(bfmix):
            print('Failing for', f)
            raise Exception("Number of basic functionals and their weights do not match.")
        for line in functionals:
            if f == line.split()[0]:
                break
        else:
            print('Failing for', f)
            raise Exception("Composite in functional_definitions.dat is missing from composite_functionals.dat")
        for bf in bfs:
            for line in basics:
                if bf == line.split()[3]:
                    if line.split()[2] == 'LIBXC' and \
                            (impl != 'LIBXC' and impl != 'BOTH'):
                        print('Failing for', f, bf)
                        raise Exception("Basic functional availibility does not match composite one.")
                    if line.split()[2] == 'XCFUN' and \
                            (impl != 'XCFUN' and impl != 'BOTH'):
                        print('Failing for', f, bf)
                        raise Exception("Basic functional availibility does not match composite one.")
                    if 'EITHER_OR' == impl and line.split()[2] != 'BOTH':
                        print('Failing for', f, bf)
                        raise Exception("Basic functional availibility does not match composite one.")
                    break
            else:
                print('Failing for', bf)
                raise Exception(
                    "Basic functional in functional_definitions.dat is missing from basic_functionals.dat")
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('Functional resolveLibXC(FUNCTIONALS functional) {\n')
    cf.write('  using BasicFunctionals::BASIC_FUNCTIONALS;\n')
    cf.write('  switch (functional) {\n')
    for d in definitions:
        split = d.split('|')
        f = split[0].split()[0]
        impl = split[0].split()[1]
        bfs = split[1].split()
        bfmix = split[2].split()
        factors = split[3].split()
        if impl in ['EITHER_OR', 'LIBXC']:
            cf.write('    case FUNCTIONALS::{}:\n'.format(f))
            cf.write('      return Functional(IMPLEMENTATIONS::LIBXC, {{{}}}, {{{}}}, {}, {}, {}, {}, {}, {});\n'.format(
                ', '.join(['BASIC_FUNCTIONALS::'+i for i in bfs]), ', '.join(bfmix), *factors))
    cf.write('    default:\n')
    cf.write('      throw SerenityError("Composite functional unknown to LibXC.");\n')
    cf.write('      break;\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('Functional resolveXCFun(FUNCTIONALS functional) {\n')
    cf.write('  using BasicFunctionals::BASIC_FUNCTIONALS;\n')
    cf.write('  switch (functional) {\n')
    for d in definitions:
        split = d.split('|')
        f = split[0].split()[0]
        impl = split[0].split()[1]
        bfs = split[1].split()
        bfmix = split[2].split()
        factors = split[3].split()
        if impl in ['EITHER_OR', 'XCFUN']:
            cf.write('    case FUNCTIONALS::{}:\n'.format(f))
            cf.write('      return Functional(IMPLEMENTATIONS::XCFUN, {{{}}}, {{{}}}, {}, {}, {}, {}, {}, {});\n'.format(
                ', '.join(['BASIC_FUNCTIONALS::'+i for i in bfs]), ', '.join(bfmix), *factors))
    cf.write('    default:\n')
    cf.write('      throw SerenityError("Composite functional unknown to XCFun.");\n')
    cf.write('      break;\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n\n')
    cf.write('} /* namespace CompositeFunctionals */\n')
    cf.write('namespace Options {\n')
    cf.write('template<>\n')
    cf.write(
        'void resolve<CompositeFunctionals::FUNCTIONALS>(std::string& value, CompositeFunctionals::FUNCTIONALS& field) {\n')
    cf.write('  static const std::map<std::string, CompositeFunctionals::FUNCTIONALS> m = {\n')
    first = True
    for i, f in enumerate(functionals):
        data = f.split()
        for alias in data[2:]:
            if not first:
                cf.write(',\n')
            first = False
            cf.write('    {{"{}", CompositeFunctionals::FUNCTIONALS::{}}}'.format(alias, data[0]))
    cf.write('\n')
    cf.write('  };\n')
    cf.write('  check(m, value, field);\n')
    cf.write('}\n\n')
    cf.write('template<>\n')
    cf.write(
        'void resolve<CompositeFunctionals::XCFUNCTIONALS>(std::string& value, CompositeFunctionals::XCFUNCTIONALS& field) {\n')
    cf.write('  static const std::map<std::string, CompositeFunctionals::XCFUNCTIONALS> m = {\n')
    first = True
    for i, f in enumerate(functionals):
        data = f.split()
        if not data[1] == 'XC' and not data[1] == 'NONE':
            continue
        for alias in data[2:]:
            if not first:
                cf.write(',\n')
            first = False
            cf.write('    {{"{}", CompositeFunctionals::XCFUNCTIONALS::{}}}'.format(alias, data[0]))
    cf.write('\n')
    cf.write('  };\n')
    cf.write('  check(m, value, field);\n')
    cf.write('}\n\n')
    cf.write('template<>\n')
    cf.write(
        'void resolve<CompositeFunctionals::KINFUNCTIONALS>(std::string& value, CompositeFunctionals::KINFUNCTIONALS& field) {\n')
    cf.write('  static const std::map<std::string, CompositeFunctionals::KINFUNCTIONALS> m = {\n')
    first = True
    for i, f in enumerate(functionals):
        data = f.split()
        if not data[1] == 'K' and not data[1] == 'NONE':
            continue
        for alias in data[2:]:
            if not first:
                cf.write(',\n')
            first = False
            cf.write('    {{"{}", CompositeFunctionals::KINFUNCTIONALS::{}}}'.format(alias, data[0]))
    cf.write('\n')
    cf.write('  };\n')
    cf.write('  check(m, value, field);\n')
    cf.write('}\n')
    cf.write('} /* namespace Options */\n\n')
    cf.write('} /* namespace Serenity */\n')

# ============= #
#   Unittests   #
# ============= #

with open('BasicFunctionals_test.cpp', 'w') as bf:
    bf.write('/**\n')
    bf.write(' * @file   BasicFunctionals_test.cpp\n')
    bf.write(' *\n')
    bf.write(' * @date   Sep 3, 2020\n')
    bf.write(' * @author Jan P. Unsleber\n')
    bf.write(' *\n')
    bf.write(' * IMPORTANT:\\n\n')
    bf.write(' * This file was automatically generated, please do not alter it.\n')
    bf.write(' * Any required changes should be made to the generating Python script\n')
    bf.write(' * which should be located close by.\n')
    bf.write(' *\n')
    bf.write(' * @copyright \\n\n')
    bf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    bf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    bf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    bf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    bf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    bf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    bf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    bf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    bf.write(' *  GNU General Public License for more details.\\n\\n\n')
    bf.write(' *  You should have received a copy of the GNU Lesser General\n')
    bf.write(' *  Public License along with Serenity.\n')
    bf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    bf.write(' */\n')
    bf.write('/* Include Class Header*/\n')
    bf.write('#include "dft/functionals/BasicFunctionals.h"\n')
    bf.write('/* Include Serenity Internal Headers */\n')
    bf.write('#include "misc/SerenityError.h"\n')
    bf.write('/* Include Std and External Headers */\n')
    bf.write('#include <gtest/gtest.h>\n')
    bf.write('\n')
    bf.write('namespace Serenity {\n')
    bf.write('namespace BasicFunctionals {\n\n')
    bf.write('#ifdef SERENITY_USE_LIBXC\n')
    bf.write('TEST(BasicFunctionalsTest, BasicFunctionalsAliases_LibXC) {\n')
    bf.write('  const unsigned int nBFuncs = {:d};\n'.format(nBFuncs))
    bf.write('  for (unsigned int i = 0; i < nBFuncs; i++){\n')
    bf.write('    auto func = static_cast<BASIC_FUNCTIONALS>(i);\n')
    bf.write('    if (func == BASIC_FUNCTIONALS::NONE || func == BASIC_FUNCTIONALS::XC_SAOP)\n')
    bf.write('      continue;\n')
    bf.write('    if (getImplementation[i] == IMPLEMENTATIONS::LIBXC || getImplementation[i] == IMPLEMENTATIONS::BOTH){\n')
    bf.write('      ASSERT_NO_THROW(getLibXCAlias(func));\n')
    bf.write('    } else {\n')
    bf.write('      ASSERT_THROW(getLibXCAlias(func), SerenityError);\n')
    bf.write('    }\n')
    bf.write('  }\n')
    bf.write('}\n')
    bf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    bf.write('#ifdef SERENITY_USE_XCFUN\n')
    bf.write('TEST(BasicFunctionalsTest, BasicFunctionalsAliases_XCFun) {\n')
    bf.write('  const unsigned int nBFuncs = {:d};\n'.format(nBFuncs))
    bf.write('  for (unsigned int i = 0; i < nBFuncs; i++){\n')
    bf.write('    auto func = static_cast<BASIC_FUNCTIONALS>(i);\n')
    bf.write('    if (func == BASIC_FUNCTIONALS::NONE || func == BASIC_FUNCTIONALS::XC_SAOP)\n')
    bf.write('      continue;\n')
    bf.write('    if (getImplementation[i] == IMPLEMENTATIONS::XCFUN || getImplementation[i] == IMPLEMENTATIONS::BOTH){\n')
    bf.write('      ASSERT_NO_THROW(getXCFunAlias(func));\n')
    bf.write('    } else {\n')
    bf.write('      ASSERT_THROW(getXCFunAlias(func), SerenityError);\n')
    bf.write('    }\n')
    bf.write('  }\n')
    bf.write('}\n')
    bf.write('#endif /* SERENITY_USE_XCFUN */\n\n')
    bf.write('TEST(BasicFunctionalsTest, BasicFunctionalsPurposes) {\n')
    bf.write('  const unsigned int nBFuncs = {:d};\n'.format(nBFuncs))
    bf.write('  for (unsigned int i = 0; i < nBFuncs; i++){\n')
    bf.write('    ASSERT_NO_THROW(getPurpose[i]);\n')
    bf.write('    int p = static_cast<int>(getPurpose[i]);\n')
    bf.write('    ASSERT_TRUE(p < 5);\n')
    bf.write('    ASSERT_TRUE(p >= 0);\n')
    bf.write('  }\n')
    bf.write('}\n\n')
    bf.write('TEST(BasicFunctionalsTest, BasicFunctionalsClass) {\n')
    bf.write('  const unsigned int nBFuncs = {:d};\n'.format(nBFuncs))
    bf.write('  for (unsigned int i = 0; i < nBFuncs; i++){\n')
    bf.write('    ASSERT_NO_THROW(getClass[i]);\n')
    bf.write('    int c = static_cast<int>(getClass[i]);\n')
    bf.write('    ASSERT_TRUE(c < 5);\n')
    bf.write('    ASSERT_TRUE(c >= 0);\n')
    bf.write('  }\n')
    bf.write('}\n\n')
    bf.write('} /* namespace BasicFunctionals */\n')
    bf.write('} /* namespace Serenity */\n\n')

with open('CompositeFunctionals_test.cpp', 'w') as cf:
    cf.write('/**\n')
    cf.write(' * @file   CompositeFunctionals_test.cpp\n')
    cf.write(' *\n')
    cf.write(' * @date   Sep 3, 2020\n')
    cf.write(' * @author Jan P. Unsleber\n')
    cf.write(' *\n')
    cf.write(' * IMPORTANT:\\n\n')
    cf.write(' * This file was automatically generated, please do not alter it.\n')
    cf.write(' * Any required changes should be made to the generating Python script\n')
    cf.write(' * which should be located close by.\n')
    cf.write(' *\n')
    cf.write(' * @copyright \\n\n')
    cf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    cf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    cf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    cf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    cf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    cf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    cf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    cf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    cf.write(' *  GNU General Public License for more details.\\n\\n\n')
    cf.write(' *  You should have received a copy of the GNU Lesser General\n')
    cf.write(' *  Public License along with Serenity.\n')
    cf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    cf.write(' */\n')
    cf.write('/* Include Class Header*/\n')
    cf.write('#include "dft/functionals/CompositeFunctionals.h"\n')
    cf.write('/* Include Serenity Internal Headers */\n')
    cf.write('#include "dft/Functional.h"\n')
    cf.write('#include "dft/functionals/BasicFunctionals.h"\n')
    cf.write('#include "misc/SerenityError.h"\n')
    cf.write('#include "settings/DFTOptions.h"\n')
    cf.write('/* Include Std and External Headers */\n')
    cf.write('#include <gtest/gtest.h>\n')
    cf.write('\n')
    cf.write('namespace Serenity {\n')
    cf.write('namespace CompositeFunctionals {\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveEnum_Functional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR){\n')
    cf.write('      ASSERT_NO_THROW(resolveFunctional(func));\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC){\n')
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('      ASSERT_NO_THROW(resolveFunctional(func));\n')
    cf.write('#else /* SERENITY_USE_LIBXC */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(func), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN){\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('      ASSERT_NO_THROW(resolveFunctional(func));\n')
    cf.write('#else /* SERENITY_USE_XCFUN */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(func), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n')
    cf.write('    }\n')
    cf.write('    else {\n')
    cf.write('      ASSERT_THROW(resolveFunctional(func), SerenityError);\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveEnum_XCFunctional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    if (getPurpose[i] != PURPOSES::EXCHANGE_CORRELATION){\n')
    cf.write('      continue;\n')
    cf.write('    }\n')
    cf.write('    auto xcfunc = static_cast<XCFUNCTIONALS>(i);\n')
    cf.write('    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR){\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC){\n')
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));\n')
    cf.write('#else /* SERENITY_USE_LIBXC */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(xcfunc), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN){\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));\n')
    cf.write('#else /* SERENITY_USE_XCFUN */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(xcfunc), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n')
    cf.write('    }\n')
    cf.write('    else {\n')
    cf.write('      ASSERT_THROW(resolveFunctional(func), SerenityError);\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveEnum_KinFunctional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    if (getPurpose[i] != PURPOSES::KINETIC){\n')
    cf.write('      continue;\n')
    cf.write('    }\n')
    cf.write('    auto kinfunc = static_cast<KINFUNCTIONALS>(i);\n')
    cf.write('    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR){\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC){\n')
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));\n')
    cf.write('#else /* SERENITY_USE_LIBXC */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(kinfunc), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n')
    cf.write('    }\n')
    cf.write('    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN){\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));\n')
    cf.write('#else /* SERENITY_USE_XCFUN */\n')
    cf.write('      ASSERT_THROW(resolveFunctional(kinfunc), SerenityError);\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n')
    cf.write('    }\n')
    cf.write('    else {\n')
    cf.write('      ASSERT_THROW(resolveFunctional(func), SerenityError);\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('#ifdef SERENITY_USE_LIBXC\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveEnum_Functional_LibXC) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    if (getImplementation[i] == IMPLEMENTATIONS::LIBXC || getImplementation[i] == IMPLEMENTATIONS::EITHER_OR){\n')
    cf.write('      ASSERT_NO_THROW(resolveLibXC(func));\n')
    cf.write('    } else {\n')
    cf.write('      ASSERT_THROW(resolveLibXC(func), SerenityError);\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('}\n')
    cf.write('#endif /* SERENITY_USE_LIBXC */\n\n')
    cf.write('#ifdef SERENITY_USE_XCFUN\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveEnum_Functional_XCFun) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    if (getImplementation[i] == IMPLEMENTATIONS::XCFUN || getImplementation[i] == IMPLEMENTATIONS::EITHER_OR){\n')
    cf.write('      ASSERT_NO_THROW(resolveXCFun(func));\n')
    cf.write('    } else {\n')
    cf.write('      ASSERT_THROW(resolveXCFun(func), SerenityError);\n')
    cf.write('    }\n')
    cf.write('  }\n')
    cf.write('}\n')
    cf.write('#endif /* SERENITY_USE_XCFUN */\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveString_Functional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    auto func = static_cast<FUNCTIONALS>(i);\n')
    cf.write('    std::string dummy;\n')
    cf.write('    FUNCTIONALS test;\n')
    cf.write('    Options::resolve<FUNCTIONALS>(dummy, func);\n')
    cf.write('    Options::resolve<FUNCTIONALS>(dummy, test);\n')
    cf.write('    ASSERT_EQ(func, test);\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveString_XCFunctional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    if (getPurpose[i] != PURPOSES::EXCHANGE_CORRELATION){\n')
    cf.write('      continue;\n')
    cf.write('    }\n')
    cf.write('    auto func = static_cast<XCFUNCTIONALS>(i);\n')
    cf.write('    std::string dummy;\n')
    cf.write('    XCFUNCTIONALS test;\n')
    cf.write('    Options::resolve<XCFUNCTIONALS>(dummy, func);\n')
    cf.write('    Options::resolve<XCFUNCTIONALS>(dummy, test);\n')
    cf.write('    ASSERT_EQ(func, test);\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('TEST(CompositeFunctionalsTest, ResolveString_KinFunctional) {\n')
    cf.write('  const unsigned int nFuncs = {:d};\n'.format(nFuncs))
    cf.write('  for (unsigned int i = 0; i < nFuncs; i++){\n')
    cf.write('    if (getPurpose[i] != PURPOSES::KINETIC){\n')
    cf.write('      continue;\n')
    cf.write('    }\n')
    cf.write('    auto func = static_cast<KINFUNCTIONALS>(i);\n')
    cf.write('    std::string dummy;\n')
    cf.write('    KINFUNCTIONALS test;\n')
    cf.write('    Options::resolve<KINFUNCTIONALS>(dummy, func);\n')
    cf.write('    Options::resolve<KINFUNCTIONALS>(dummy, test);\n')
    cf.write('    ASSERT_EQ(func, test);\n')
    cf.write('  }\n')
    cf.write('}\n\n')
    cf.write('} /* namespace CompositeFunctionals */\n')
    cf.write('} /* namespace Serenity */\n\n')

# =================== #
#   Python Bindings   #
# =================== #

with open('BasicFunctionals_python.cpp', 'w') as bf:
    bf.write('/**\n')
    bf.write(' * @file   BasicFunctionals_python.cpp\n')
    bf.write(' *\n')
    bf.write(' * @date   Feb 18, 2025\n')
    bf.write(' * @author Anton Rikus\n')
    bf.write(' *\n')
    bf.write(' * IMPORTANT:\\n\n')
    bf.write(' * This file was automatically generated, please do not alter it.\n')
    bf.write(' * Any required changes should be made to the generating Python script\n')
    bf.write(' * which should be located close by.\n')
    bf.write(' *\n')
    bf.write(' * @copyright \\n\n')
    bf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    bf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    bf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    bf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    bf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    bf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    bf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    bf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    bf.write(' *  GNU General Public License for more details.\\n\\n\n')
    bf.write(' *  You should have received a copy of the GNU Lesser General\n')
    bf.write(' *  Public License along with Serenity.\n')
    bf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    bf.write(' */\n')
    bf.write('/* Include Serenity Internal Headers */\n')
    bf.write('#include "dft/functionals/BasicFunctionals.h"\n')
    bf.write('/* Include Std and External Headers */\n')
    bf.write('#include <pybind11/pybind11.h>\n\n')
    bf.write('using namespace Serenity;\n')
    bf.write('namespace py = pybind11;\n\n')
    bf.write('void export_BasicFunctionals(py::module& m) {\n\n')
    bf.write('  py::enum_<BasicFunctionals::BASIC_FUNCTIONALS>(m, "BASIC_FUNCTIONALS")\n')
    for i, f in enumerate(basics):
        data = f.split()
        bf.write('      .value("{}", BasicFunctionals::BASIC_FUNCTIONALS::{})\n'.format(data[3], data[3]))
    bf.write('      .export_values();\n')
    bf.write('}\n')

with open('CompositeFunctionals_python.cpp', 'w') as cf:
    cf.write('/**\n')
    cf.write(' * @file   CompositeFunctionals_python.cpp\n')
    cf.write(' *\n')
    cf.write(' * @date   Sep 3, 2020\n')
    cf.write(' * @author Jan P. Unsleber\n')
    cf.write(' *\n')
    cf.write(' * IMPORTANT:\\n\n')
    cf.write(' * This file was automatically generated, please do not alter it.\n')
    cf.write(' * Any required changes should be made to the generating Python script\n')
    cf.write(' * which should be located close by.\n')
    cf.write(' *\n')
    cf.write(' * @copyright \\n\n')
    cf.write(' *  This file is part of the program Serenity.\\n\\n\n')
    cf.write(' *  Serenity is free software: you can redistribute it and/or modify\n')
    cf.write(' *  it under the terms of the GNU Lesser General Public License as\n')
    cf.write(' *  published by the Free Software Foundation, either version 3 of\n')
    cf.write(' *  the License, or (at your option) any later version.\\n\\n\n')
    cf.write(' *  Serenity is distributed in the hope that it will be useful,\n')
    cf.write(' *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    cf.write(' *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    cf.write(' *  GNU General Public License for more details.\\n\\n\n')
    cf.write(' *  You should have received a copy of the GNU Lesser General\n')
    cf.write(' *  Public License along with Serenity.\n')
    cf.write(' *  If not, see <http://www.gnu.org/licenses/>.\\n\n')
    cf.write(' */\n')
    cf.write('/* Include Serenity Internal Headers */\n')
    cf.write('#include "dft/functionals/CompositeFunctionals.h"\n')
    cf.write('/* Include Std and External Headers */\n')
    cf.write('#include <pybind11/pybind11.h>\n\n')
    cf.write('using namespace Serenity;\n')
    cf.write('namespace py = pybind11;\n\n')
    cf.write('void export_CompositeFunctionals(py::module& spy) {\n\n')
    cf.write('  py::enum_<CompositeFunctionals::FUNCTIONALS>(spy, "FUNCTIONALS")\n')
    for i, f in enumerate(functionals):
        data = f.split()
        cf.write('      .value("{}", CompositeFunctionals::FUNCTIONALS::{})\n'.format(data[2], data[0]))
    cf.write('      .export_values();\n')
    cf.write('  py::enum_<CompositeFunctionals::XCFUNCTIONALS>(spy, "XCFUNCTIONALS")\n')
    for i, f in enumerate(functionals):
        data = f.split()
        if data[1] == 'XC' or data[1] == 'NONE':
            cf.write('      .value("{}", CompositeFunctionals::XCFUNCTIONALS::{})\n'.format(data[2], data[0]))
    cf.write('      .export_values();\n')
    cf.write('  py::enum_<CompositeFunctionals::KINFUNCTIONALS>(spy, "KINFUNCTIONALS")\n')
    for i, f in enumerate(functionals):
        data = f.split()
        if data[1] == 'K' or data[1] == 'NONE':
            cf.write('      .value("{}", CompositeFunctionals::KINFUNCTIONALS::{})\n'.format(data[2], data[0]))
    cf.write('      .export_values();\n')
    cf.write('}\n')
