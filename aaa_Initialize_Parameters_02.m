function aaa_Initialize_Parameters_02
% Get screen resolution
scrn_res_pxl = get(0, "ScreenSize"); % Screen Resolution Pixe
% Create a window for placing the items for user input
wh = scrn_res_pxl(4)/1.4; % Window Height
ww = 600; % Window Width
win_left = (scrn_res_pxl(3) - ww)/2; % Window Left Location
win_bottom = (scrn_res_pxl(4) - wh)/2; % Window Bottom Location
f = uifigure("Position", [win_left win_bottom ww wh]); % Creat Window

g = uigridlayout(f);
n_row = 20;
for i = 1 : n_row
    rh_cell(i) = {'fit'};
end 
clear i
g.RowHeight = rh_cell;
g.ColumnWidth = {'fit', '1x'};

clr_ws = uicheckbox(g);
clr_ws.Text = "Clear workspace?"; 
clr_ws.Value = 1;
clr_ws.Layout.Row = 1;
clr_ws.Layout.Column = 1;

save_mat = uicheckbox(g);
save_mat.Text = "Save a mat file with all the essential variables?"; 
save_mat.Value = 1;
save_mat.Layout.Row = 2;
save_mat.Layout.Column = 1;

manual_refine = uicheckbox(g);
manual_refine.Text = "Manually refine the ROI selection?"; 
manual_refine.Value = 0;
manual_refine.Layout.Row = 3;
manual_refine.Layout.Column = 1;

crop_image_panel = uipanel(g);
crop_image_panel.Title = "Cropping";
crop_image_panel.Layout.Row = 4;
crop_image_panel.Layout.Column = [1,2];
crop_image_grid = uigridlayout(crop_image_panel);
crop_image_grid.RowHeight = {'fit', 'fit'};
crop_image_grid.ColumnWidth = {'fit', '1x'};

crop_image = uicheckbox(crop_image_grid);
crop_image.Text = "Crop the image stack?"; 
crop_image.Value = 0;
crop_image.Layout.Row = 1;
crop_image.Layout.Column = 1;

crop_coordinate_text = uilabel(crop_image_grid);
crop_coordinate_text.Text = " X1 | X2 | Y1 | Y2"; 
crop_coordinate_text.Layout.Row = 2;
crop_coordinate_text.Layout.Column = 2;
%->
crop_coordinate = uieditfield(crop_image_grid);
crop_coordinate.Value = "100, 450, 300, 450"; 
crop_coordinate.Layout.Row = 2;
crop_coordinate.Layout.Column = 1;

do_scn_phs_crcn = uicheckbox(g);
do_scn_phs_crcn.Text = "Do line scan phase correction?"; 
do_scn_phs_crcn.Value = 1;
do_scn_phs_crcn.Layout.Row = 6;
do_scn_phs_crcn.Layout.Column = 1;

do_mtn_crcn = uicheckbox(g);
do_mtn_crcn.Text = "Do motion correction?"; 
do_mtn_crcn.Value = 1;
do_mtn_crcn.Layout.Row = 7;
do_mtn_crcn.Layout.Column = 1;

mk_mv = uicheckbox(g);
mk_mv.Text = "Save the processed image stack?"; 
mk_mv.Value = 0;
mk_mv.Layout.Row = 8;
mk_mv.Layout.Column = 1;

n_avrg_text = uilabel(g);
n_avrg_text.Text = "Number of average"; 
n_avrg_text.Layout.Row = 9;
n_avrg_text.Layout.Column = 2;
%->
n_avrg = uieditfield(g);
n_avrg.Value = "6"; 
n_avrg.Layout.Row = 9;
n_avrg.Layout.Column = 1;

brk_panel = uipanel(g);
brk_panel.Title = "Breaking the stack";
brk_panel.Layout.Row = 10;
brk_panel.Layout.Column = [1,2];
brk_grid = uigridlayout(brk_panel);
brk_grid.RowHeight = {'fit', 'fit', 'fit'};
brk_grid.ColumnWidth = {'fit', '1x'};

brk_stck = uicheckbox(brk_grid);
brk_stck.Text = "Break the stack to smaller stacks?"; 
brk_stck.Value = 0;
brk_stck.Layout.Row = 1;
brk_stck.Layout.Column = 1;
%->
n_sgmt_text = uilabel(brk_grid);
n_sgmt_text.Text = "Number of sub-stacks (segments)"; 
n_sgmt_text.Layout.Row = 2;
n_sgmt_text.Layout.Column = 2;
%->
n_sgmt = uieditfield(brk_grid);
n_sgmt.Value = "3"; 
n_sgmt.Layout.Row = 2;
n_sgmt.Layout.Column = 1;
%->
trgt_sgmt_text = uilabel(brk_grid);
trgt_sgmt_text.Text = "Which sub-stack to process"; 
trgt_sgmt_text.Layout.Row = 3;
trgt_sgmt_text.Layout.Column = 2;
%->
trgt_sgmt = uieditfield(brk_grid);
trgt_sgmt.Value = "1"; 
trgt_sgmt.Layout.Row = 3;
trgt_sgmt.Layout.Column = 1;

gssn_xyz_text = uilabel(g);
gssn_xyz_text.Text = "X | Y | Z | 3D Guassian smoothing parameters"; 
gssn_xyz_text.Layout.Row = 13;
gssn_xyz_text.Layout.Column = 2;
%->
gssn_xyz = uieditfield(g);
gssn_xyz.Value = "1, 0.9, 1"; 
gssn_xyz.Layout.Row = 13;
gssn_xyz.Layout.Column = 1;

hstm_stch_tshd_text = uilabel(g);
hstm_stch_tshd_text.Text = "Histogram streching threshhold"; 
hstm_stch_tshd_text.Layout.Row = 14;
hstm_stch_tshd_text.Layout.Column = 2;
%->
hstm_stch_tshd = uieditfield(g);
hstm_stch_tshd.Value = "2000"; 
hstm_stch_tshd.Layout.Row = 14;
hstm_stch_tshd.Layout.Column = 1;

gssn_hstm_stch_tshd_text = uilabel(g);
gssn_hstm_stch_tshd_text.Text = "Histogram streching threshhold after smoothing"; 
gssn_hstm_stch_tshd_text.Layout.Row = 15;
gssn_hstm_stch_tshd_text.Layout.Column = 2;
%->
gssn_hstm_stch_tshd = uieditfield(g);
gssn_hstm_stch_tshd.Value = "5000"; 
gssn_hstm_stch_tshd.Layout.Row = 15;
gssn_hstm_stch_tshd.Layout.Column = 1;

%--CaImAn Setting Pannel------------
caiman_panel = uipanel(g);
caiman_panel.Title = "CaImAn Settings";
caiman_panel.Layout.Row = 16;
caiman_panel.Layout.Column = [1,2];
caiman_grid = uigridlayout(caiman_panel);
caiman_grid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
caiman_grid.ColumnWidth = {'fit', '1x'};
%-->
phs_crcn_frm_text = uilabel(caiman_grid);
phs_crcn_frm_text.Text = "Phase correction number of frames "; 
phs_crcn_frm_text.Layout.Row = 1;
phs_crcn_frm_text.Layout.Column = 2;
%->
phs_crcn_frm = uieditfield(caiman_grid);
phs_crcn_frm.Value = "200"; 
phs_crcn_frm.Layout.Row = 1;
phs_crcn_frm.Layout.Column = 1;
%-->
nrn_half_text = uilabel(caiman_grid);
nrn_half_text.Text = "Neuron half size in pixels"; 
nrn_half_text.Layout.Row = 2;
nrn_half_text.Layout.Column = 2;
%->
nrn_half = uieditfield(caiman_grid);
nrn_half.Value = "5"; 
nrn_half.Layout.Row = 2;
nrn_half.Layout.Column = 1;
%-->
n_nrn_text = uilabel(caiman_grid);
n_nrn_text.Text = "Number of neurons (components)"; 
n_nrn_text.Layout.Row = 3;
n_nrn_text.Layout.Column = 2;
%->
n_nrn = uieditfield(caiman_grid);
n_nrn.Value = "40"; 
n_nrn.Layout.Row = 3;
n_nrn.Layout.Column = 1;
%-->
ca_prfl_text = uilabel(caiman_grid);
ca_prfl_text.Text = "Calcium transient profile; 1: decay; 2: rise & decay"; 
ca_prfl_text.Layout.Row = 4;
ca_prfl_text.Layout.Column = 2;
%->
ca_prfl = uieditfield(caiman_grid);
ca_prfl.Value = "2"; 
ca_prfl.Layout.Row = 4;
ca_prfl.Layout.Column = 1;
%-->
bg_cmpt_text = uilabel(caiman_grid);
bg_cmpt_text.Text = "Number of background components"; 
bg_cmpt_text.Layout.Row = 5;
bg_cmpt_text.Layout.Column = 2;
%->
bg_cmpt = uieditfield(caiman_grid);
bg_cmpt.Value = "1"; 
bg_cmpt.Layout.Row = 5;
bg_cmpt.Layout.Column = 1;


% ---------------------------
% Exit Callback 
exitReturn = uibutton(g);
exitReturn.Layout.Row = n_row;
exitReturn.Layout.Column = 1;
exitReturn.Text = "Save to workspace";
exitReturn.ButtonPushedFcn = @exitReturn_Callback;
uiwait(f);

    function exitReturn_Callback(hObject, eventdata)
        disp("Bye Bye")
        assignin("base", "clr_ws", clr_ws.Value);
        assignin("base", "save_mat", save_mat.Value);
        assignin("base", "save_mat", save_mat.Value);
        assignin("base", "manual_refine", manual_refine.Value);
        assignin("base", "crop_image", crop_image.Value);
        assignin("base", "crop_coordinate", str2num(crop_coordinate.Value));
        assignin("base", "do_mtn_crcn", do_mtn_crcn.Value);
        assignin("base", "do_scn_phs_crcn", do_scn_phs_crcn.Value);
        assignin("base", "mk_mv", mk_mv.Value);
        assignin("base", "n_avrg", str2num(n_avrg.Value));
        assignin("base", "brk_stck", brk_stck.Value);
        assignin("base", "n_sgmt", str2num(n_sgmt.Value));
        assignin("base", "trgt_sgmt", str2num(trgt_sgmt.Value));
        assignin("base", "gssn_xyz", str2num(gssn_xyz.Value));
        assignin("base", "hstm_stch_tshd", str2num(hstm_stch_tshd.Value));
        assignin("base", "gssn_hstm_stch_tshd", str2num(gssn_hstm_stch_tshd.Value));
        assignin("base", "phs_crcn_frm", str2num(phs_crcn_frm.Value));
        assignin("base", "n_nrn", str2num(n_nrn.Value));
        assignin("base", "nrn_half", str2num(nrn_half.Value));
        assignin("base", "ca_prfl", str2num(ca_prfl.Value));
        assignin("base", "bg_cmpt", str2num(bg_cmpt.Value));
        close(f);
    end

end 