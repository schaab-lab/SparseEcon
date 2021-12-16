clear
close all
addpath(genpath("../lib/"))
delete(gcp('nocreate')); parpool;

diary ./output/unit_tests.log
diary on

tic
test_setup_grid;
test_gen_H;
test_gen_BH;
test_parents;
test_gen_FD;
test_keep_dims;
test_adapt_grid;
toc

diary off
close all