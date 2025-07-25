%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-20 22:08:17
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:28:40
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% Test script for MatCell class

% Create MatCell objects from cell arrays
cell_input1 = {rand(2, 2), rand(2, 2)};
cell_input2 = {rand(2, 2), rand(2, 2)};

mc1 = MatCell(cell_input1);
mc2 = MatCell(cell_input2);

fprintf('display MatCell object 1:\n');
disp(mc1);

% Access a block using subsref
block1 = mc1{1};
assert(isequal(block1, cell_input1{1}), 'Error: subsref failed');

% Assign a block using subsasgn
new_block = rand(2, 2);
mc3 = mc1;
mc3{1} = new_block;
assert(isequal(mc3{1}, new_block), 'Error: subsasgn failed');


% Perform arithmetic operations on MatCell objects
mc_sum = mc1 + mc2;
mc_diff = mc1 - mc2;
mc_neg = -mc1;
mc_mult = mc1 * mc2;
mc_dotmult = mc1 .* mc2;
mc_div = mc1 ./ mc2;

% Check arithmetic operations results
assert(isequal(mc_sum{1}, cell_input1{1} + cell_input2{1}), 'Error: plus operation failed');
assert(isequal(mc_diff{1}, cell_input1{1} - cell_input2{1}), 'Error: minus operation failed');
assert(isequal(mc_neg{1}, -cell_input1{1}), 'Error: uminus operation failed');
assert(isequal(mc_mult{1}, cell_input1{1} * cell_input2{1}) && isequal(mc_mult{2}, cell_input1{2} * cell_input2{2}), 'Error: mtimes operation failed');
assert(isequal(mc_dotmult{1}, cell_input1{1} .* cell_input2{1}), 'Error: times operation failed');
assert(isequal(mc_div{1}, cell_input1{1} ./ cell_input2{1}), 'Error: rdivide operation failed');


cell_input1 = {rand(2, 3), rand(2, 3)};
mc1 = MatCell(cell_input1);


% Reduce sum of a MatCell object
mc_reduce_sum = MatCell.reduce_sum(mc1);
assert(isequal(mc_reduce_sum, cell_input1{1} + cell_input1{2}), 'Error: reduce_sum operation failed');

% Vertical concatenation of a MatCell object
mc_vert_concat = MatCell.vert_concat(mc1);
assert(isequal(mc_vert_concat, [cell_input1{1}; cell_input1{2}]), 'Error: vert_concat operation failed');

% Horizontal concatenation of a MatCell object
mc_hori_concat = MatCell.hori_concat(mc1);
assert(isequal(mc_hori_concat, [cell_input1{1}, cell_input1{2}]), 'Error: hori_concat operation failed');


% Create a matrix A
A = rand(4, 6);

% Define the block column sizes in an array cone_size
cone_size = [2, 4];


% Split the matrix A into blocks using the vert_split method
mc_split = MatCell.vert_split(A', cone_size);


% Split the matrix A into blocks using the vert_split method
mc_split = MatCell.hori_split(A, cone_size);


% Display success message
disp('All tests passed successfully');