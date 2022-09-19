0 my_help called with 'list ruby -f'
1 - ruby help.
2 - [[https://syntaxdb.com/ref/ruby][Syntax_DB]] ruby syntac 
3 - http://rubular.com regular expression
4 - [[https://github.com/DragonRuby/command_line][command_line]] idea for command line test
5 - [[https://github.com/peterc/testrocket][testrocket]] idea for simple test
6 -----
7 [0;32;49mfsm_finite_state_machine[0m
8 #+begin_src ruby
9     require 'scanf'
10 
11     file = 'relax_calc.o22733'
12 
13     TRANS = {
14       wait_next: {
15         '* fix ' => [:searching, :xy_data],
16         :default => [:wait_next, :idle]
17       },
18       searching: {
19         '   1 F' => [:wait_next, :z_data],
20         :default => [:searching, :idle]
21       }
22     }
23     state = :wait_next
24     data = []
25     all_data = []
26     File.readlines(file).each do |line|
27       state, action = TRANS[state][line[0..5]] || TRANS[state][:default]
28       case action
29       when :idle
30         data = line.scanf("* fix calc kpoints:50, in_plane:%f, vertical:%f")
31       when :z_data
32         data << line.scanf("   1 F= %f E0= %f  d E =%f")[0]
33         all_data << data.flatten
34       end
35     end
36 
37     pp all_data
38   # https://qiita.com/daddygongon/private/67c8d40c276fdec472c7
39 #+end_src
40 
