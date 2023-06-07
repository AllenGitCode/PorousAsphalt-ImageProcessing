%% Example to show useage of "detectConnectedPoreCells"
% In this example, we apply the detect-connected-pores algorithm on a
% sample of porous media. If you have access to your own porous media
% dataset, you can replace the example sample with your own sample. Details
% of the detect-connected-pores algorithm are presented in Allen et al.,
% 2023. "Development of a MATLAB-based code for quantification of effective
% void space in porous pavement." Presented at the 64th International
% Conference of Scandinavian Simulation Society, SIMS 2023, Västerås,
% Sweden, September 26-27, 2023.

addpath ../utils/
addpath ../data/

%% Loading and visualization of a dataset
% There are a few possibilities here: 1 - upload a few images, and do
% filtering, binarization, visualization on the images. 2 - make a
% synthetic sample using, for example, MATLAB's built-in rand() function. 3
% - call a function that gives a pre-coded dataset. 4 - load a .mat file
% containing a pre-processed dataset. We chose to carry out our example
% using option 4. That is, we have already applied masking, cropping,
% filtering, and binarization of hundreds of .tiff images to generate the
% dataset. We simply load the results here.
load dataset_sampleA.mat
volumeViewer(BW_matrix_stack, 'ScaleFactors', [1 1 1])


%% Plot porosity profile
% quantify porosity of each 2D image found in the stack
total_cells_in_roi = sum(sum(mask)); % total cells in ROI
for i = 1:length(sliceNums)
    BWpore = BW_matrix_stack(:,:,i);
    total_pore_cells_in_roi = sum(sum(BWpore));
    poro(i) = total_pore_cells_in_roi / total_cells_in_roi;
end
figure, hold on
plot(poro,corresponding_elevation_of_slice,'-','DisplayName','total')
xlabel('cross-sectional porosity')
ylabel('elevation (cm)')
legend('show')
xlim([0 0.3])
box


%% Pass over to the detectConnectedPoreCells algorithm
% get connected cells, visualize

% optional: reduce resolution of dataset to speed up processing time
vertRes = 50;
sliceNums_coarsened = sliceNums(1):vertRes:sliceNums(end);
corresponding_elevation_of_slice_coarsened = corresponding_elevation_of_slice(1:vertRes:end);
BW_matrix_stack_coarsened = BW_matrix_stack(:,:,1:vertRes:end);
BW_colArray_stack_coarsened = reshape(BW_matrix_stack_coarsened,[],1);

% then pass to algorithm
t = tic; disp(['processing effective pore space...'])
[connected_matrix_stack, ~, countTOTAL, countTOP, countBOTTOM] = detectConnectedPoreCells( BW_colArray_stack_coarsened, m, n );
elaspedTime = toc(t); disp(['Elasped time to process effective pore space is ' num2str(elaspedTime/60) ' minutes,'])

% then visualize
volumeViewer(connected_matrix_stack, 'ScaleFactors', [1 1 vertRes])


%% Plot connected porosity profile
% compared to total porosity profile
total_cells_in_roi = sum(sum(mask)); % total cells in ROI
for i = 1:length(sliceNums_coarsened)
    BWpore = connected_matrix_stack(:,:,i);
    total_conn_pore_cells_in_roi = sum(sum(BWpore));
    poro_connected(i) = total_conn_pore_cells_in_roi / total_cells_in_roi;
end
figure, hold on
plot(poro,corresponding_elevation_of_slice,'-','DisplayName','total')
plot(poro_connected,corresponding_elevation_of_slice_coarsened,'--','DisplayName','connected')
xlabel('cross-sectional porosity')
ylabel('elevation (cm)')
legend('show')
xlim([0 0.3])
box


%% Algorithm reporting
% report iteration number, CPU run time
table(elaspedTime/60, countTOP, countBOTTOM, countTOTAL, ...
    'VariableNames',{'CPU time (mins)','# while-loop its. for TOP-connected','# while-loop its. for BOTTOM-connected','total # while-loop its.'})



%%
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