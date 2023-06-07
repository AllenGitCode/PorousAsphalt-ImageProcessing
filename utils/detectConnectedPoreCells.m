function [connected_matrix_stack, connected_colArray, countTOTAL, countTOP, countBOTTOM] = detectConnectedPoreCells( BW_colArray_stack, m_crop, n_crop )
% Detect the pore cells that are connected from top-to-bottom. Top is the
% m_crop*n_crop BW image, where black (logical true) are pore cells and
% white (logical false) are non-pore cells (including any outside ROI
% cells).
%
% SYNOPSIS:
%   connected_matrix_stack = detectConnectedPoreCells( BW_colArray_stack, m_crop, n_crop )
%   [connected_matrix_stack, connected_colArray, ...
%       countTOTAL, countTOP, countBOTTOM] = detectConnectedPoreCells( BW_colArray_stack, m_crop, n_crop )
%
% PARAMETERS:
%   BW_colArray_stack   -   A 1D column-array containing elements of 1 and 0.
%                           Elements of 0 are solid or outside roi.
%                           Elements of 1 are pore cells.
%
%   m_crop  - Number of rows in 2D image
%
%   n_crop  - Number of columns in 2D image
%
% RETURNS:
%   connected_matrix_stack - ...
%   
%   connected_colArray - ...
%   
%   countTOTAL - ...
%   
%   countTOP - ...
%   
%   countBOTTOM - ...
%
% NOTE:
%   Function detectConnectedPoreCells assumes that
%   length(BW_colArray_stack) >= m_crop*n_crop
%
% EXAMPLES:
%   % 
%   m = 100; n = 100; numSlices = 5;
%   BW_matrix_stack = rand(m,n,numSlices);
%   BW_colArray_stack = reshape(BW_matrix_stack,[m,n,numSlices]);
%   detectConnectedPoreCells(BW_colArray_stack, m, n)
%
% SEE ALSO:
%   

%{
Copyright 2023 Rebecca Allen

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the “Software”),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
%}

%% Number of 2D slices that make up a 3D stack or volume
length_of_sliceNums = size(BW_colArray_stack,1) / (m_crop*n_crop);

%% Get cell numbers (their index) of all 3D model pore cells. We use int32
% and not uint32 to avoid error when using ismember later on.
poreCells = int32(find(BW_colArray_stack));
topPoreCells = int32(find(BW_colArray_stack(1:m_crop*n_crop)));
bottomPoreCells = int32(find(BW_colArray_stack((end-(m_crop*n_crop)+1):end)));
bottomPoreCells = bottomPoreCells + int32((m_crop*n_crop*length_of_sliceNums)-(m_crop*n_crop)+1);

%% Create a surrounding cells map, with the following structure:
%
%   - col 1 is main pore cell number
%   - col 2 is label of main pore cell number
%   - col 3 is R (right) surrounding cell number
%   - col 4 is L (left) surrounding cell number
%   - col 5 is F (front) surrounding cell number
%   - col 6 is B (back) surrounding cell number
%   - col 7 is U (up) surrounding cell number
%   - col 8 is D (down) surrounding cell number
%
% We note that a surrounding cell number might not be a pore. We also note
% that is no surrounding cell in a certain direction exists, then this
% element is assigned as outside roi, i.e., a void/pore space
%
poreCellsMap = int32(poreCells);
%
% main pore cell label, initialized with a 2 (different from 0 and 1)
poreCellsMap(:,2) = int32(2); 
%
% surrounding cell labels
poreCellsMap(:,3) = int32(poreCellsMap(:,1) + 1); % R surrounding cell number is main pore cell + 1, except at bdry
poreCellsMap(:,4) = int32(poreCellsMap(:,1) - 1); % L surrounding cell number is main pore cell - 1, except at bdry
poreCellsMap(:,5) = int32(poreCellsMap(:,1) + m_crop); % F surrounding cell number is main pore cell + m, except at bdry
poreCellsMap(:,6) = int32(poreCellsMap(:,1) - m_crop); % B surrounding cell number is main pore cell - m, except at bdry
poreCellsMap(:,7) = int32(poreCellsMap(:,1) + (m_crop*n_crop)); % U surrounding cell number is main pore cell + m*n, except at bdry
poreCellsMap(:,8) = int32(poreCellsMap(:,1) - (m_crop*n_crop)); % D surrounding cell number is main pore cell - m*n, except at bdry


%% Some of the surrounding cells are actually outside of the VOI (volume of
% interest), and we do not want to count them as possible connected pore
% cells. Thus we replace any surrounding cell number that is outside the
% VOI with a Nan, with the following three steps:
% 
% First, we get the indices of the boundary cells:
bdryCells = [];
for i = 1:length_of_sliceNums
    % South bdry cells
    bdryCells(1:m_crop,i) = (1:m_crop) + (i-1)*m_crop*n_crop;
    % North bdry cells
    bdryCells((m_crop+1):(2*m_crop),i) = ((n_crop*m_crop-m_crop+1):(n_crop*m_crop)) + (i-1)*m_crop*n_crop;
    % West bdry cells
    for j = 1:n_crop
        bdryCells(2*m_crop+j,i) = (j-1)*m_crop+1 + (i-1)*m_crop*n_crop;
    end
    % East bdry cells
    for j = 1:n_crop
        bdryCells(2*m_crop+n_crop+j,i) = (j-1)*m_crop+m_crop + (i-1)*m_crop*n_crop;
    end
end
% reshape into colArray and take only unique bdry cells
bdryCells = reshape(bdryCells,[],1);
bdryCells = unique(bdryCells);
%
% Second, we replace any elements of poreCellsMap(:,3) to poreCellsMap(:,8)
% that contain values equal to the bdryCells values with a NaN (or
% appropriate bdry flag). This means the actual surrounding cell to a main
% pore cell was outside the VOI bdry, however it was wrongly labeled
% with a +1 shifted periodic bdry cell and we want to replace it with a NaN
logicalIndices = ismember(poreCellsMap(:,3), bdryCells);
poreCellsMap(logicalIndices,3) = nan;
logicalIndices = ismember(poreCellsMap(:,4), bdryCells);
poreCellsMap(logicalIndices,4) = nan;
logicalIndices = ismember(poreCellsMap(:,5), bdryCells);
poreCellsMap(logicalIndices,5) = nan;
logicalIndices = ismember(poreCellsMap(:,6), bdryCells);
poreCellsMap(logicalIndices,6) = nan;
logicalIndices = ismember(poreCellsMap(:,7), bdryCells);
poreCellsMap(logicalIndices,7) = nan;
logicalIndices = ismember(poreCellsMap(:,8), bdryCells);
poreCellsMap(logicalIndices,8) = nan;
%
% Third, we replace any elements of poreCellsMap(:,8) and poreCellsMap(:,7)
% that contain values less than 1 or values greater than the max domain
% cell number, respectively.
logicalIndices = poreCellsMap(:,8) < 1;
poreCellsMap( logicalIndices, 8) = nan;
logicalIndices = poreCellsMap(:,7) > m_crop*n_crop*length_of_sliceNums;
poreCellsMap( logicalIndices, 7) = nan;
%
poreCellsMap_init = poreCellsMap; % for later use


%% Now that a map of surrounding cells to main pore cells is complete, we
% replace col-2 with correct "label". That is, all TOP pore cells get
% labeled as a top-connected pore cell. And all pore cells that have at
% least one surrounding cell that is labeled as a top-connected pore cell
% get labeled as a top-connected pore cell. This procedure of updating
% labels should continue until there are no "newly labeled" top-connected
% pore cells.
%
disp(['   -detecting connected pore cells from TOP'])
%
% label TOP pore cells with a 1 (TOP-connected pore cell)
indexes = int32(find(ismember(poreCellsMap(:,1),topPoreCells)));
poreCellsMap(indexes, 2) = 1;
%
c = 1;
countTOP = 0;
while ~isempty(find(c))
    countTOP = countTOP + 1;
    % get indexes of all pore cells currently labeled as TOP-connected pore cells
    indexes = int32(find(poreCellsMap(:,2)==1));
    % get these pore cell numbers
    topConnectedPoreCellNumbers = poreCellsMap(indexes,1);
    % look for these numbers in any of the surrounding cells, get index, and
    % replace poreCellsMap(:,2) with appropriate label
    c = ismember(poreCellsMap(:,3:8), topConnectedPoreCellNumbers);
    % c will eventually be empty if surrounding cells that are top-connected
    % get replaced with other value
    indexes_3 = find(c(:,3-2));
    indexes_4 = find(c(:,4-2));
    indexes_5 = find(c(:,5-2));
    indexes_6 = find(c(:,6-2));
    indexes_7 = find(c(:,7-2));
    indexes_8 = find(c(:,8-2));
    poreCellsMap(indexes_3,2) = 1;
    poreCellsMap(indexes_4,2) = 1;
    poreCellsMap(indexes_5,2) = 1;
    poreCellsMap(indexes_6,2) = 1;
    poreCellsMap(indexes_7,2) = 1;
    poreCellsMap(indexes_8,2) = 1;
    % also replace poreCellsMap(:,3:8) with other value so element is not
    % detected again in following loop iteration
    poreCellsMap(indexes_3,3) = nan;
    poreCellsMap(indexes_4,4) = nan;
    poreCellsMap(indexes_5,5) = nan;
    poreCellsMap(indexes_6,6) = nan;
    poreCellsMap(indexes_7,7) = nan;
    poreCellsMap(indexes_8,8) = nan;
end
%
TOPconnectedPoreCellsIndex = int32(find(poreCellsMap(:,2)==1));
TOPconnectedPoreCells = poreCellsMap(TOPconnectedPoreCellsIndex);
clear TOPconnectedPoreCellsIndex


%% Then, we need to do a scan from bottom to detect bottom connected cells
disp(['   -detecting connected pore cells from BOTTOM'])
M = m_crop*n_crop*length_of_sliceNums;
%
% get fresh poreCellsMap
poreCellsMap = poreCellsMap_init; clear poreCellsMap_init
%
% label BOTTOM pore cells with a 1 (BOTTOM-connected pore cell)
indexes = int32(find(ismember(poreCellsMap(:,1),bottomPoreCells)));
poreCellsMap(indexes, 2) = 1;
%
% repeat while loop
c = 1;
countBOTTOM = 0;
while ~isempty(find(c))
    countBOTTOM = countBOTTOM + 1;
    % get indexes of all pore cells currently labeled as TOP-connected pore cells
    indexes = int32(find(poreCellsMap(:,2)==1));
    % get these pore cell numbers
    bottomConnectedPoreCellNumbers = poreCellsMap(indexes,1);
    % look for these numbers in any of the surrounding cells, get index, and
    % replace poreCellsMap(:,2) with appropriate label
    c = ismember(poreCellsMap(:,3:8), bottomConnectedPoreCellNumbers);
    % c will eventually be empty if surrounding cells that are bottom-connected
    % get replaced with other value
    indexes_3 = find(c(:,3-2)); 
    indexes_4 = find(c(:,4-2));
    indexes_5 = find(c(:,5-2));
    indexes_6 = find(c(:,6-2));
    indexes_7 = find(c(:,7-2));
    indexes_8 = find(c(:,8-2));
    poreCellsMap(indexes_3,2) = 1;
    poreCellsMap(indexes_4,2) = 1;
    poreCellsMap(indexes_5,2) = 1;
    poreCellsMap(indexes_6,2) = 1;
    poreCellsMap(indexes_7,2) = 1;
    poreCellsMap(indexes_8,2) = 1;
    % also replace poreCellsMap(:,3:8) with other value so element is not
    % detected again in following loop iteration
    poreCellsMap(indexes_3,3) = nan;
    poreCellsMap(indexes_4,4) = nan;
    poreCellsMap(indexes_5,5) = nan;
    poreCellsMap(indexes_6,6) = nan;
    poreCellsMap(indexes_7,7) = nan;
    poreCellsMap(indexes_8,8) = nan;
end
%
BOTTOMconnectedPoreCellsIndex = int32(find(poreCellsMap(:,2)==1));
BOTTOMconnectedPoreCells = poreCellsMap(BOTTOMconnectedPoreCellsIndex);
clear BOTTOMconnectedPoreCellsIndex


%% Then, we find all members between top and bottom connected cells, to get
% the top-to-bottom connected pore cells
disp(['   -merging connected pore cells from TOP and BOTTOM'])
%
connectedPoreCellsIndex = int32(find(ismember(TOPconnectedPoreCells, BOTTOMconnectedPoreCells)));
clear BOTTOMconnectedPoreCells
connectedPoreCells = TOPconnectedPoreCells(connectedPoreCellsIndex);
clear connectedPoreCellsIndex
clear TOPconnectedPoresCells
%
connected_colArray = int8(zeros(M, 1));
connected_colArray(connectedPoreCells) = 1;
connected_matrix_stack = reshape(connected_colArray, [m_crop,n_crop,length_of_sliceNums]);
%
if max(max(max(connected_matrix_stack))) == 0
    error(['Error: \nsample with selected slice numbers and roi radius...' ...
        ' \ndoes not result in a top-to-bottom connected pore space. ...' ...
        ' \nTry using a different range of slice numbers and radius.'])
end


%% Lastly, we get isolated pore cells (optional)
disp(['   -extracting isolated pore cells'])
isolated_colArray = int8(BW_colArray_stack); % 1 = pore space, 0 = solid or outside roi
isolated_colArray(connectedPoreCells) = 0;
isolated_matrix_stack = reshape(isolated_colArray, [m_crop,n_crop,length_of_sliceNums]);
clear isolated_colArray
clear BW_colArray_stack
clear connectedPoreCells


%% For reporting purposes
countTOTAL = countTOP + countBOTTOM;


end