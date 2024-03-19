Contribute to Serenity
===========================

- [Contribute to Serenity](#contribute-to-serenity)
    - [Contributor License Agreement](#contributor-license-agreement)
    - [Who Can Use My Unpublished Algorithms](#who-can-use-my-unpublished-algorithms)
    - [Issue Tracker and Feature Requests](#issue-tracker)
    - [Merge Requests](#merge-requests)
        - [Contribution Acceptance Criteria](#contribution-acceptance-criteria)
    - [Code and Style Guidelines](#code-and-style-guidelines)
        - [Object Names](#object-names)
        - [Namespaces](#namespaces)
        - [Include Statements](#include-statements)
        - [Whitespaces and Line Breaks](#whitespaces-and-line-breaks)
        - [Object Comments](#object-comments)
        - [Algorithm Comments](#algorithm-comments)
        - [Language](#language)
        - [Copies, References and Pointers](#copies,-references-and-pointers)
    - [Code of conduct](#code-of-conduct)


## Contributor License Agreement
By submitting code as an individual you agree to the usage and sharing of
that code under the GNU L-GPLv3 as provided in the LICENSE file.
If you are a developer, in addition to this please read the section:
[Who Can Use My Unpublished Algorithms](#who-can-use-my-unpublished-algorithms)

## Access
While this code is intended to be open source, we are aware that peoples' academic careers might be 
intertwined with the code they add to this project.
The intention here is obviously for everyone to publish their own algorithms/methods before handing
them to other people.
For this reason, the development git repository is not open to everyone.  
General users of the code are kindly asked to use the release versions hosted on
[GitHub](https://github.com/qcserenity/serenity).
Every developer will be invited and given access to the development repository and thus be able to clone 
the entire repository with all branches.

## Who Can Use My Unpublished Algorithms

This section pertains only to code that is still exclusively present in the closed development repository.

The general rule is that all changes merged into the `master` are free to be used by any developer.
They are considered ready for release and ready for production tests.
The usage of algorithms under development in un-merged branches for publication should be off 
limits unless agreed upon by the developing person(s).  
Infringement of this last point can and will result in a permanent removal from the project's development team.  
Additionally, the citation of key publications in reference with Serenity is mandatory.

## Issue Tracker and Feature Requests
Search the issue tracker for similar entries before submitting your own.
Use apropriate labels please.
Keep feature proposals as small and simple as possible, complex ones might be edited to make them small and simple.
Obviously, some features need chunks of new code - please use tasks to break down bigger feature proposals.

## Merge Requests
In general, all improvements to the code are welcome, however, in order to get your
merge request accepted, please read the following sections carefully.

### Contribution Acceptance Criteria

Your merge will be accepted if you:
1. Include proper tests and make all tests pass (unless it contains a test
   exposing a bug in existing code). In general, every new class or function 
   should have sensible corresponding unit tests.
2. If you suspect a failing test is unrelated to your contribution, please ask about it.
3. Make sure your changes can merge without problems (otherwise, merge the
   current `master` branch back into your branch).
4. Do not break any existing functionality.
5. Keep the Serenity code base clean and well structured
   There is a simple bash script that can format the code using `clang-format`.
   It is located under `dev/scripts/formatting.sh`.
6. Conform to the [Code and Style Guidelines](#code-and-style-guidelines)
7. Have reasonable documentation. If we can't get what you did by looking at the
   documentation, "You shall not pass!".
8. Give yourself credit, add a Changelog entry if it is more than a 'one-line' bug fix.

## Code and Style Guidelines
  
There are many rules and guidelines out there that describe how to write clean code,
going through them all would be beyond the scope of this document.
In addition, a lot comes down to personal style, however, on top of the few small
guidelines given in the following sections please always try to do the following.
Generally think about how you write your code. Ask yourself:
   - Will other people understand your code?
   - Will you understand your code yourself after not looking into it for some time?
   - Is the code readable and sensible?
   - May names you used be easily misinterpreted?

To improve you can always add comments and re-choose/elongate names you chose.

Oh yeah, and please, try to be consistent with other parts of the project.

### Object Names
Class names: always start with upper case. New words inside the name are indicated by
upper case characters. The name of a class is equal to the name of the files
defining it and giving the method implementations.
  - Example: ExampleClassName
    (which is defined in 'ExampleClassName.h' with method implementations
    in 'ExampleClassName.cpp')

Public variable names: start with a lower case letter, followed by mixed case
  - Example: exampleVariableName

Private variable names: are indicated by a preceding underscore.
  - Example: _privateMemberVariable

Global/Static constant names: have upper case letters only. Words are seperated by an underscore.
  - Examples: GLOBAL_CONSTANT, PI

Function names: are named like normal variables, i.e. mixed case starting with a lower
case letter.
  - Example: exampleFunction()

### Namespaces
Everything inside this project is defined within the namespace 'Serenity'.
Subnamespaces are allowed. Names are to be chosen not to conflict with
namespace 'std'.


### Include Statements
Avoid includes in header files.
Before you put an include statement into a header file, try to use a forward declaration 
It is most often sufficient to have the include in the .cpp file.
Includes should be ordered, most of the time it is done like this:
The first include is reserved for the header file corresponding to the
current source file. Then all other Serenity-internal includes.
Finally includes for all external projects and the standard libraries follow.  
Inside each block includes are ordered alphabetically. Includes which are possibly temporary and
for miscellaneous things like Timing or Print routines may be placed in a separate block.

### Whitespaces and Line Breaks
Line Length should be limited at 120 characters.
We use two white spaces for indentation. An exception is the namespace Serenity;
no indentation is done here.
NO tabs, please!


### Object Comments
Please use comments which can be interpreted by doxygen for everything which is
not too detailed. As a rough rule: Make all comments (for public variables and
members, as well as general class and file comments) in header files interpretable
by doxygen (start with \/\*\*) and none (except for file comments) in source files.
The \@ symbol is used for doxygen commands.

Every file has a little comment in the top specifying at least the date of file
creation and your full real name.

### Algorithm Comments
On top of the class, function, file, etc. documentation, there is always the
good old documentation inside the functions.
Please document your algorithms, noone likes numbers falling from the sky.
Arguably good variable names do a lot of the work, but when in doubt noone
will complain about too much documentation.

### Language
...is English. This is not only valid for comments, but also names you define
should make sense to an English speaking person. 

### Copies, References and Pointers

There could be a lot said here, but really it boils down to the folowing few guidelines
for this package and its data:
  - Keep things as 'const' as possible
  - Try to copy as little as possible
  - Use smart pointers, avoid raw pointers at almost any cost

## Code of conduct

As contributors and maintainers of this project, we pledge to respect all
people who contribute through reporting issues, posting feature requests,
updating documentation, submitting pull requests or patches, and other
activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual
language or imagery, derogatory comments or personal attacks, trolling, public
or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct. Project maintainers who do not
follow the Code of Conduct may be removed from the project team.

This code of conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

This Code of Conduct is adapted from the Contributor Covenant, version 1.1.0,
available at [http://contributor-covenant.org/version/1/1/0/](http://contributor-covenant.org/version/1/1/0/).
  
