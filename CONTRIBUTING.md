# Contributing to Project_2025

## Styleguides

In order to ease the development and mantainance of the code, please follow the
following styleguides. This will make sure that the code is readable, not only
for you but for the rest of developers.

In case something is not specified below, take as a reference the code that is
currently present in this repository.

- Comment your code and name your variables in english.

- Use descriptive variable names. (E.g.: `positions` or `pos` are acceptable
  variable names for a positions array. `x`, `p`, `pstns`, `posi`, `asdfghjk` or
  `potato` are not).

- Use proper code indentation. Do loops, if statements, subroutines, functions,
  etc. must always be indented.

  Example of acceptable code indentation:

```fortran
module io
    implicit none

contains
    ! Read data from a specified .xyz file.
    subroutine read_xyz(input_file, molecule_name, num_points, atom, x_axis, y_axis, z_axis)
        implicit none

        character(*), intent(in) :: input_file
        character(*), intent(out) :: molecule_name
        character(2), allocatable, intent(out) :: atom(:)
        real, allocatable, intent(out) :: x_axis(:), y_axis(:), z_axis(:)
        integer, intent(out) :: num_points
        integer :: i

        open(10, file = input_file, status = 'old')

        write(*, '(3a)') 'Reading values from ', trim(adjustl(input_file)), '...'
        write(*, *)

        read(10, *) num_points
        read(10, '(a)') molecule_name

        allocate(atom(num_points))
        allocate(x_axis(num_points))
        allocate(y_axis(num_points))
        allocate(z_axis(num_points))

        do i = 1, num_points
            read(10, *) atom(i), x_axis(i), y_axis(i), z_axis(i)
        end do

        close(10)
    end subroutine read_xyz

end module io
```

  Examples of bad code indentation:

```fortran
module abc
implicit none

contains
subroutine xyz(a, b, c)
implicit none

...

do i = 1, x
x = x + i
end do
end subroutine xyz
end module abc
```

Example of unacceptable code indentation:

```fortran
module abc
                implicit none
        contains
    subroutine xyz
implicit none
...
```

- Really long lines of code are bad for code readability and code clarity,
  worsening the overall code quality and making maintenance more difficult.

  Avoid lines extending over 90 characters of length and __never make lines
  longer than 132 characters__. Lines longer than 132 characters might cause
  compilation errors and/or unexpected code behaviour

- Comment your code but not over-comment it.

  Code comments at the beginning of a function or subroutine or in order to
  explain or clarify certain non-trivial parts of the code are helpful and
  necessary. Filling the code with unnecessary comments can (and will) difficult
  code comprehension for everyone.

- Avoid leaving commented chunks of code. If you do some tests in your code that
  you end up disabling by commenting them out and are sure they won't be used
  again, remove them instead of leaving them commented. Big chunks of commented
  code, again, worsen code readability.

- Avoid the use of capital letters in your code and __never write in all caps__.
  Fortran doesn't distinguish between caps or non-caps characters. So all
  non-caps variables like `positions_output_file` are preferred over
  `PositionsOutputFile` or `POSITIONSOUTPUTFILE`, as they are better
  readability-wise.

- In case you need to add a module to the program, add it in a separate file.
  There must only be one module per file.

## Commit Messages

- Please, separate your code changes and updates into several different commits
  in order to be able to better track your changes. It is prefereable to have
  several small commits where each of them is related to a small specific change
  in the code that one large commit containing all your changes. Spliting your
  changes in different commits helps with the tracking of changes and can help
  us see when a bug was first introduced and solve the issues it may cause
  faster.

- Commits must have a short and descriptive title that can give an overall idea
  of the changes you are submitting.

- It is also a good practice to add a short description about what the commit
  contains in order to provide more context on the submitted changes (although
  it is not necessary for commits that only add files or do trivial changes).

Example of good commits (first line is the title and the rest of the body is the
description):

```
Add Makefile
```

```
Add run_temperature_loop.sh

Add a shell script to perform iterative runs of the full program code
for different temperatures.
```

```
Add skip_steps input parameter

Change the skip_lines variable (previously hard coded to 1) to a more
appropiate name: header_lines.

As the user might want to skip a series of Monte Carlo steps when
computing binning (E.g.: in order to only take into account the
equilibrated system), the skip_steps input variable is added. This
variable is defined by the user in 'input_parameters.in'.
```

Example of bad commits and commit messages:

```
Modified abc.f90 due to an error with the xyz. Modified defg.f90 to correct the way the hijkl array was being accessed and modified the conditional for the...
```

Notice how the commit text is only in the title and the commit contains several
unrelated changes. Different unrelated changes must be separated into different
commits in order to be able to better track changes.

```
Changes
```

This commit provides no information about the changes.
