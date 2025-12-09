#include "SetParameterWindow.h"
#include "FramedPanel.h"
#include "RolloutPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

SetParameterWindow::SetParameterWindow() : BaseWindow(600, 800, "Set parameter")
{
    // Setup for drawing パネルを作成します。
    create_vertical_space(15);

    // Setup for drawing パネルのサイズです。
    const int setup_width = 550;
    const int setup_height = 400;

    RolloutPanel* setup_panel = new RolloutPanel(
        window_width, setup_width, setup_height, "Setup for drawing", true);
    setup_panel->begin();

    // Setup for drawing パネルの 1 行目の項目群です。
    const int setup_row1_width = setup_width - 30;
    const int setup_row1_height = 90;
    Fl_Pack* setup_row1 = create_row(setup_row1_width, setup_row1_height, 10);

    // Setup for drawing パネルの 1 行目の 1 列目です。
    const int setup_row1_column1_width = 140;
    Fl_Pack* setup_row1_column1 = create_column(setup_row1_column1_width, setup_row1_height);

    // 原子のヒット機能のオンオフを行うチェック ボックスを追加します。2024/11 CJS
    create_check_button(setup_row1_column1_width, "Enable atom click", &enable_atom_click, 1);

    create_check_button(setup_row1_column1_width, "Ortho", &ortho_state, 0);

    setup_row1_column1->end();

    // Setup for drawing パネルの 1 行目の 2 列目です。
    const int setup_row1_column2_width = 100;
    Fl_Pack* setup_row1_column2 = create_column(setup_row1_column2_width, setup_row1_height);

    create_check_button(setup_row1_column1_width, "Show axis", &show_axis);
    create_check_button(setup_row1_column1_width, "Show cell", &show_cell);

    setup_row1_column2->end();

    // Setup for drawing パネルの 1 行目の 3 列目です。
    const int setup_row1_column3_width = 220;
    Fl_Pack* setup_row1_column3 = create_column(setup_row1_column3_width, setup_row1_height);

    create_int_spinner(setup_row1_column3_width, 80, "Sphere radius:", &radius, 1, 20, 5);
    create_int_spinner(setup_row1_column3_width, 80, "Sphere segment:", &segments, 0, 200);
    create_float_spinner(setup_row1_column3_width, 80, "Replica range:", &replica_range, 0.0, 0.5, 0.01);

    setup_row1_column3->end();
    setup_row1->end();

    // Setup for drawing パネルの 2 行目の項目群です。
    const int setup_panel_row2_width = 500;
    const int setup_panel_row2_height = 190;
    Fl_Pack* setup_panel_row2 = create_row(
        setup_panel_row2_width, setup_panel_row2_height, 10);

    // Bond パネルを作成します。
    int bond_panel_width = 300;
    int bond_panel_height = 160;

    Fl_Group* bond_frame = create_framed_group_xy(
        0, 0, bond_panel_width, bond_panel_height, "Bond");

    Fl_Pack* bond_pack = create_panel_pack_xy(
        0, 30, bond_panel_width, bond_panel_height, false);

    // 1 列目です。
    // 列の幅です。
    const int bond_column1_width = 100;
    Fl_Pack* bond_column1 = create_column(bond_column1_width, bond_panel_height);

    create_vertical_space(10);
    create_check_button(bond_column1_width, "Draw bond", &draw_bond, CB_BOND);
    create_check_button(bond_column1_width, "Bond-PBC", &draw_bond_pbc, CB_BOND);
    create_radio_button(bond_column1_width, "Line", &bond_display_type, 0, CB_BOND);
    create_radio_button(bond_column1_width, "Cylinder", &bond_display_type, 1, CB_BOND);
    create_dropdown_list(80, 80, "", &bond_display_color, color_list, MAXCOLOR, CB_BOND);

    bond_column1->end();

    // 2 列目です。
    const int bond_column2_width = 200;
    Fl_Pack* bond_column2 = create_column(bond_column2_width, bond_panel_height);

    create_vertical_space(10);
    create_float_spinner(bond_column2_width - 20, 90, "Bond length:", &bondlength, 0.0, 10.0, 0.1, CB_BOND);
    create_int_spinner(bond_column2_width - 20, 90, "Bond thick:", &bondthickness, 0, 20, CB_BOND);

    bond_column2->end();
    bond_pack->end();
    bond_frame->end();

    // Force/Load パネルです。
    int force_load_panel_width = 180;
    int force_load_panel_height = 120;
    Fl_Group* force_load_frame = create_framed_group_xy(
        -10, 0, force_load_panel_width, force_load_panel_height, "Force/Load");

    Fl_Pack* force_load_panel = create_panel_pack_xy(
        0, 30, force_load_panel_width, force_load_panel_height, true);

    create_check_button(force_load_panel_width, "Draw Force", &draw_force, 10);
    create_check_button(force_load_panel_width, "Draw Load", &draw_load, 10);
    create_float_input(force_load_panel_width - 20, 60, "F-arw thick:", &vscl_force);
    create_check_button(force_load_panel_width, "Draw Auxiliary", &draw_aux, 10);

    force_load_panel->end();
    force_load_frame->end();
    setup_panel_row2->end();

    create_vertical_space(10);

    // Setup for drawing パネルの 3 行目の項目です。
    Fl_Pack* setup_row3 = create_row(setup_width, 30);
    create_horizontal_space(setup_width / 2 - 60);
    create_int_spinner(150, 60, "Redraw itvl:", &mdspeed, 1, 50);
    setup_row3->end();
    setup_panel->end();

    // Atom color パネルを作成します
    // Atom color パネルのサイズです。
    const int atom_color_panel_width = 500;
    const int atom_color_panel_height = 230;
    RolloutPanel* atom_color_panel = new RolloutPanel(
        window_width, atom_color_panel_width, atom_color_panel_height, "Atom color", false);
    atom_color_panel->spacing(10);
    atom_color_panel->begin();

    // Atom color パネルの 1 列目です。
    const int atom_color_column1_width = 120;
    Fl_Pack* atom_color_column1 = create_column(
        atom_color_column1_width, atom_color_panel_height);
    create_vertical_space(10);
    // ラジオ ボタン グループを作成します。
    create_radio_button(120, "Species", &color_mode, 0, CB_COLOR_MODE);
    create_radio_button(120, "Energy", &color_mode, 1, CB_COLOR_MODE);
    create_radio_button(120, "CSP", &color_mode, 2, CB_COLOR_MODE);
    create_radio_button(120, "Mises", &color_mode, 3, CB_COLOR_MODE);

    // 入力フィールドを作成します。
    create_float_input(120, 60, "min:", &color_mode_vmin, CB_COLOR_MODE);
    create_float_input(120, 60, "max:", &color_mode_vmax, CB_COLOR_MODE);

    // チェック ボックスを作成します。
    create_check_button(120, "Autoscale", &color_mode_auto, CB_COLOR_MODE);
    atom_color_column1->end();

    create_horizontal_space(20);

    // Atom color パネルの 2 列目です。
    const int atom_color_column2_width = 200;
    Fl_Pack* atom_color_column2 = create_column(
        atom_color_column2_width, atom_color_panel_height);
    create_vertical_space(10);
    create_int_spinner(atom_color_column2_width, 60, "Marked atom 0:", &marked_atom[0], 0, 10000, CB_MARKED_ATOM);
    create_int_spinner(atom_color_column2_width, 60, "Marked atom 1:", &marked_atom[1], 0, 10000, CB_MARKED_ATOM);
    create_int_spinner(atom_color_column2_width, 60, "Marked atom 2:", &marked_atom[2], 0, 10000, CB_MARKED_ATOM);
    atom_color_column2->end();

    // Atom color パネルの 3 列目です。
    const int atom_color_column3_width = 100;
    Fl_Pack* atom_color_column3 = create_column(
        atom_color_column3_width, atom_color_panel_height);
    create_vertical_space(10);
    // カラー選択のドロップダウン リストです。
    create_dropdown_list(
        60, 60, "", &marked_atom_color[0], color_list, MAXCOLOR, CB_MARKED_ATOM);
    create_dropdown_list(
        60, 60, "", &marked_atom_color[1], color_list, MAXCOLOR, CB_MARKED_ATOM);
    create_dropdown_list(
        60, 60, "", &marked_atom_color[2], color_list, MAXCOLOR, CB_MARKED_ATOM);
    atom_color_column3->end();

    atom_color_panel->end();

    // Cutoff/Bookkeep パネルを作成します。
    // Cutoff/Bookkeep パネルのサイズです。
    const int cutoff_bookkeep_panel_width = 550;
    const int cutoff_bookkeep_panel_height = 110;

    RolloutPanel* cutoff_bookkeep_panel = new RolloutPanel(
        window_width, cutoff_bookkeep_panel_width, cutoff_bookkeep_panel_height, "Cutoff/Bookeep", false);
    cutoff_bookkeep_panel->spacing(10);
    cutoff_bookkeep_panel->begin();

    // Cutoff/Bookkeep パネルの 1 列目です。
    const int cutoff_bookkeep_column1_width = 140;

    Fl_Pack* cutoff_bookkeep_column1 = create_column(
        cutoff_bookkeep_column1_width, cutoff_bookkeep_panel_height);

    create_vertical_space(10);
    create_float_input(cutoff_bookkeep_column1_width, 60, "Cutoff (A):", &rcut_f);

    cutoff_bookkeep_column1->end();

    // Cutoff/Bookkeep パネルの 2 列目です。
    const int cutoff_bookkeep_column2_width = 180;
    Fl_Pack* cutoff_bookkeep_column2 = create_column(
        cutoff_bookkeep_column2_width, cutoff_bookkeep_panel_height);

    create_vertical_space(10);
    create_float_input(cutoff_bookkeep_column2_width, 60, "B-keep margin (A):", &frcmar);
    create_float_input(cutoff_bookkeep_column2_width, 60, "B-keep cutoff (A):", &frc_f);

    cutoff_bookkeep_column2->end();

    // Cutoff/Bookkeep パネルの 3 列目です。
    const int cutoff_bookkeep_column3_width = 140;
    Fl_Pack* cutoff_bookkeep_column3 = create_column(
        cutoff_bookkeep_column3_width, cutoff_bookkeep_panel_height);

    create_vertical_space(10);
    create_int_spinner(cutoff_bookkeep_column3_width, 60, "B-keep itvl:", &book.nbk);

    cutoff_bookkeep_column3->end();
    cutoff_bookkeep_panel->end();

    // Set parameters ウィンドウの 4 行目を作成します。
    int setparam_panel_row4_width = 500;
    int setparam_panel_row4_height = 240;

    Fl_Pack* setpanel_row4 = create_row(setparam_panel_row4_width, setparam_panel_row4_height, 10);

    // Set parameters ウィンドウの 4 行目 1 列目です。
    const int setpanel_row4_column1_width = 180;
    Fl_Pack* setpanel_row4_column1 = create_column(setpanel_row4_column1_width, setparam_panel_row4_height);

    // Relax algo パネルです。
    const int relax_algo_panel_width = setpanel_row4_column1_width;
    const int relax_algo_panel_height = 150;

    Fl_Group* relax_algo_frame = create_framed_group_xy(
        0, 0, relax_algo_panel_width, relax_algo_panel_height, "Relax algo");

    Fl_Pack* relax_algo_pack = create_panel_pack_xy(
        0, 30, relax_algo_panel_width, relax_algo_panel_height, true);

    // ラジオ ボタンを作成します。
    create_radio_button(relax_algo_panel_width, "GLOC", &relax_algo, 0, CB_REL);
    create_radio_button(relax_algo_panel_width, "FIRE", &relax_algo, 1, CB_REL);
    create_radio_button(relax_algo_panel_width, "CG", &relax_algo, 2, CB_REL);
    create_float_input(relax_algo_panel_width - 40, 40, "Damper on V:", &relax_damper_value);

    relax_algo_pack->end();
    relax_algo_frame->end();

    create_float_input(relax_algo_panel_width, 60, "NH mass x 10^:", &nhmass_scale);

    setpanel_row4_column1->end();

    // 列 1、2 の間のスペースを調整します。
    create_horizontal_space(10);

    // Set parameters ウィンドウの 4 行目 2 列目です。
    const int setpanel_row4_column2_width = 150;
    Fl_Pack* setpanel_row4_column2 = create_column(setpanel_row4_column2_width, setparam_panel_row4_height);
    create_vertical_space(30);
    create_int_input(setpanel_row4_column2_width, 60, "CONF Wr itvl:", &confwrint);
    create_check_button(setpanel_row4_column2_width, "Auto capture", &autocap);
    setpanel_row4_column2->end();

    // 列 2、3 の間のスペースを調整します。
    create_horizontal_space(10);

    // Set parameters ウィンドウの 4 行目 3 列目です。
    int setpanel_row4_column3_width = 160;

    Fl_Pack* setpanel_row4_column3 = create_column(setpanel_row4_column3_width, setparam_panel_row4_height);
    create_vertical_space(30);
    create_check_button(setpanel_row4_column3_width, "Follow recipe", &irecipe);
    create_check_button(setpanel_row4_column3_width, "Stop MD by Fmax", &itolfor);
    create_float_input(setpanel_row4_column3_width, 50, "Fmax tlrc (eV/A):", &tolfor);
    create_check_button(setpanel_row4_column3_width, "^-- Reheat after stop", &ireheat);
    create_check_button(setpanel_row4_column3_width, "Stop MD by step", &itolstep);
    create_int_input(setpanel_row4_column3_width, 50, "# of steps:", &tolstep);
    create_check_button(setpanel_row4_column3_width, "Stop MD by stress", &itolstress);
    create_float_input(setpanel_row4_column3_width, 50, "Stress tlrc (MPa):", &tolstress);

    setpanel_row4_column3->end();

    setpanel_row4->end();

    // Reheat param パネルを作成します。
    // Reheat param パネルのサイズです。
    const int reheat_param_panel_width = 400;
    const int reheat_param_panel_height = 120;

    RolloutPanel* reheat_param_panel = new RolloutPanel(
        window_width, reheat_param_panel_width, reheat_param_panel_height, "Reheat param", false);
    reheat_param_panel->begin();

    // 列の幅です。
    const int reheat_column_width = 160;

    // Reheat param パネルの 1 列目です。
    Fl_Pack* reheat_column1 = create_column(
        reheat_column_width, reheat_param_panel_height);
    create_vertical_space(10);
    create_float_input(reheat_column_width, 80, "Temp:", &reheat_temp);
    create_float_input(reheat_column_width, 80, "dt (fs):", &reheat_dtm);
    reheat_column1->end();

    create_horizontal_space(20);

    // Reheat param パネルの 2 列目です。
    Fl_Pack* reheat_column2 = create_column(
        reheat_column_width, reheat_param_panel_height);
    create_vertical_space(10);
    create_int_input(reheat_column_width, 40, "# of steps:", &reheat_step);
    create_int_input(reheat_column_width, 40, "# of reheats to go:", &reheat_count);
    reheat_column2->end();

    reheat_param_panel->end();

    // Fire param パネルを作成します。
    const int fire_param_panel_width = 550;
    const int fire_param_panel_height = 120;

    RolloutPanel* fire_param_panel = new RolloutPanel(
        window_width, fire_param_panel_width, fire_param_panel_height, "Fire param", false);
    fire_param_panel->begin();

    // 列の幅です。
    const int fire_param_column_width = fire_param_panel_width / 3 - 10;

    create_horizontal_space(10);

    // Fire param パネルの 1 列目です。
    Fl_Pack* fire_param_column1 = create_column(fire_param_column_width, fire_param_panel_height);
    create_vertical_space(10);
    create_float_input(fire_param_column_width, 50, "alpha init:", &fp_alph_ini, CB_FIRE);
    create_int_input(fire_param_column_width, 50, "min accel steps:", &fp_nfmin, CB_FIRE);
    fire_param_column1->end();

    create_horizontal_space(10);

    // Fire param パネルの 2 列目です。
    Fl_Pack* fire_param_column2 = create_column(fire_param_column_width - 30, fire_param_panel_height);
    create_vertical_space(10);
    create_float_input(fire_param_column_width - 30, 50, "dt inc (>1):", &fp_ffinc, CB_FIRE);
    create_float_input(fire_param_column_width - 30, 50, "dt inc (<1):", &fp_ffdec, CB_FIRE);
    fire_param_column2->end();

    create_horizontal_space(10);

    // Fire param パネルの 3 列目です。
    Fl_Pack* fire_param_column3 = create_column(fire_param_column_width - 30, fire_param_panel_height);
    create_vertical_space(10);
    create_float_input(fire_param_column_width - 30, 50, "al rate (<1):", &fp_ffalph, CB_FIRE);
    create_float_input(fire_param_column_width - 30, 50, "dt max (fs):", &fp_ffdtmax, CB_FIRE);
    fire_param_column3->end();

    fire_param_panel->end();

    // Full rlx (sta) パネルを作成します。
    const int full_rlx_panel_width = 550;
    const int full_rlx_panel_height = 100;

    RolloutPanel* full_rlx_panel = new RolloutPanel(
        window_width, full_rlx_panel_width, full_rlx_panel_height, "Full rlx (sta)", true);
    full_rlx_panel->begin();

    Fl_Pack* full_rlx_row1 = create_row(full_rlx_panel_width, 30, 20);
    create_float_input(200, 60, "Fmax for cell adjust:", &cell_relax_tolfor);
    create_int_input(250, 60, "# of cell adjusts each time:", &cell_relax_rep);
    full_rlx_row1->end();
    full_rlx_panel->end();

    // Deformation settings パネルを作成します。
    const int deformation_setting_panel_width = 550;
    const int deformation_setting_panel_height = 565;

    RolloutPanel* deformation_setting_panel = new RolloutPanel(
        window_width, deformation_setting_panel_width, deformation_setting_panel_height, "Deformation settings", true);
    deformation_setting_panel->begin();

    // Grab grouped atoms パネル作成 (YU2025)
    FramedPanel* grabgroup_panel = new FramedPanel(
        deformation_setting_panel_width, deformation_setting_panel_width - 60, 255, "Grab/push grouped atoms", true);
    grabgroup_panel->begin();
    Fl_Pack* grabgroup_column1 = create_column(deformation_setting_panel_width - 60, 205);
    //Fl_Pack* grabgroup_column1 = create_column(150,180);
    Fl_Pack* grabgroup_row1 = create_row(deformation_setting_panel_width/20, 25);
    create_check_button(150, "Grab 1st group", &ifgrab1);
    create_int_spinner(150, 80, "Group #:", &grab1num, 0, MAXGROUP);
    grabgroup_row1->end();
    Fl_Pack* grabgroup_row2 = create_row(deformation_setting_panel_width/20, 25);
    create_float_input(180, 60, "dx,dy,dz (A/ns)", &grabdxdt1);
    create_float_input( 80, 60,                "", &grabdydt1);
    create_float_input( 80, 60,                "", &grabdzdt1);
    grabgroup_row2->end();
    Fl_Pack* grabgroup_row3 = create_row(deformation_setting_panel_width/20, 25);
    create_check_button(150, "Push 1st group", &ifpush1);
    create_int_spinner(150, 80, "Group #:", &push1num, 0, MAXGROUP);
    grabgroup_row3->end();
    Fl_Pack* grabgroup_row4 = create_row(deformation_setting_panel_width/20, 25);
    create_float_input(180, 60, "fx,fy,fz (nN)", &pushfx1);
    create_float_input( 80, 60,              "", &pushfy1);
    create_float_input( 80, 60,              "", &pushfz1);
    grabgroup_row4->end();
    //grabgroup_column1->end();
    //create_horizontal_space(50);
    //Fl_Pack* grabgroup_column2 = create_column(150,180);
    Fl_Pack* grabgroup_row5 = create_row(deformation_setting_panel_width/20, 25);
    create_check_button(150, "Grab 2nd group", &ifgrab2);
    create_int_spinner(150, 80, "Group #:", &grab2num, 0, MAXGROUP);
    grabgroup_row5->end();
    Fl_Pack* grabgroup_row6 = create_row(deformation_setting_panel_width/20, 25);
    create_float_input(180, 60, "dx,dy,dz (A/ns)", &grabdxdt2);
    create_float_input( 80, 60,                "", &grabdydt2);
    create_float_input( 80, 60,                "", &grabdzdt2);
    grabgroup_row6->end();
    Fl_Pack* grabgroup_row7 = create_row(deformation_setting_panel_width/20, 25);
    create_check_button(150, "Push 2nd group", &ifpush2);
    create_int_spinner(150, 80, "Group #:", &push2num, 0, MAXGROUP);
    grabgroup_row7->end();
    Fl_Pack* grabgroup_row8 = create_row(deformation_setting_panel_width/20, 25);
    create_float_input(180, 60, "fx,fy,fz (nN)", &pushfx2);
    create_float_input( 80, 60,              "", &pushfy2);
    create_float_input( 80, 60,              "", &pushfz2);
    grabgroup_row8->end();
    //grabgroup_column2->end();
    grabgroup_column1->end();
    grabgroup_panel->end();

    // Stretch/Shrink パネルを作成します。
    FramedPanel* stretch_panel = new FramedPanel(
        deformation_setting_panel_width, deformation_setting_panel_width - 60, 120, "Stretch/Shrink", true);
    stretch_panel->begin();

    // 行を追加していきます。
    Fl_Pack* stretch_row1 = create_row(deformation_setting_panel_width - 60, 30, 20);
    create_float_input(120, 60, "ex(/ps):", &dexdt);
    create_float_input(120, 60, "ey(/ps):", &deydt);
    create_float_input(120, 60, "ez(/ps):", &dezdt);
    stretch_row1->end();

    create_vertical_space(10);

    Fl_Pack* stretch_row2 = create_row(deformation_setting_panel_width - 60, 30, 20);
    create_check_button(120, "Repeat Lz", &repeat_lz);
    create_float_input(120, 60, "Lz(min):", &repeat_lz_min);
    create_float_input(120, 60, "Lz(max):", &repeat_lz_max);
    stretch_row2->end();

    stretch_panel->end();

    // Auxiliary パネルを作成します。
    FramedPanel* auxiliary_panel = new FramedPanel(
        deformation_setting_panel_width, deformation_setting_panel_width - 60, 90, "Auxiliary", true);
    auxiliary_panel->begin();

    // 水平方向のパックを作成します。
    Fl_Pack* auxiliary_row = create_row(deformation_setting_panel_width - 60, 30, 20);

    create_check_button(120, "yz plane punch", &yzplane_punch);
    create_float_input(120, 60, "d (A):", &yzplane_punch_d);
    create_float_input(160, 60, "d rate (A/ps):", &yzplane_punch_dd);

    auxiliary_row->end();
    auxiliary_panel->end();

    deformation_setting_panel->end();

    // Special settings for CNT パネルを作成します。
    const int cnt_setting_panel_width = 550;
    const int cnt_setting_panel_height = 400;

    RolloutPanel* cnt_setting_panel = new RolloutPanel(
        window_width, cnt_setting_panel_width, cnt_setting_panel_height, "Special settings for CNT", true);
    cnt_setting_panel->begin();

    // Special setting for CNT パネルの 1 行目です。
    Fl_Pack* cnt_settings_row1 = create_panel_pack(
        window_width, cnt_setting_panel_width, cnt_setting_panel_height - 100, false);
    cnt_settings_row1->spacing(10);

    // 列を追加します。
    const int cnt_settings_column_width = cnt_setting_panel_width / 2 - 10;
    const int cnt_settings_column_height = cnt_setting_panel_height - 70;
    Fl_Pack* cnt_settings_column1 = create_column(cnt_settings_column_width, cnt_settings_column_height);
    create_vertical_space(10);
    create_check_button(180, "Corrugation mode", &mode_cnt_corrugation);
    create_check_button(180, "Show CNT wall", &show_cnt_wall);
    create_check_button(180, "Show CNT n-vec", &show_cnt_wallv);
    create_check_button(180, "Show CNT ring", &show_cnt_ring);
    create_float_input(150, 60, "n-vec length:", &vscl);

    // Loading type パネルです。
    const int loading_type_panel_width = cnt_settings_column_width - 20;
    const int loading_type_panel_height = 120;

    Fl_Group* loading_type_frame = create_framed_group_xy(
        0, 0, loading_type_panel_width - 50, loading_type_panel_height, "Loading type");
    Fl_Pack* loading_type_pack = create_panel_pack_xy(
        10, 30, loading_type_panel_width - 50, loading_type_panel_height, true);

    // ラジオ ボタンを作成します。
    create_radio_button(loading_type_panel_width - 20, "Hydro (Wall)", &cnt_load_algo, 1);
    create_radio_button(loading_type_panel_width - 20, "Ring-force", &cnt_load_algo, 2);
    create_radio_button(loading_type_panel_width - 20, "Ring-rigit", &cnt_load_algo, 3);
    loading_type_pack->end();
    loading_type_frame->end();

    cnt_settings_column1->end();

    Fl_Pack* cnt_settings_column2 = create_column(cnt_setting_panel_width / 2 - 50, cnt_setting_panel_height - 70);
    create_vertical_space(10);

    // Wall パネルです。
    const int wall_panel_width = cnt_settings_column_width;
    const int wall_panel_height = 100;

    Fl_Group* wall_frame = create_framed_group_xy(
        0, 0, wall_panel_width - 50, wall_panel_height, "Wall");
    Fl_Pack* wall_pack = create_panel_pack_xy(
        10, 30, wall_panel_width - 50, wall_panel_height, true);

    // スピナーを作成します。
    create_int_spinner(200, 100, "CNT wall #:", &show_cnt_wall_num, 0, 10000, 5);
    create_float_spinner(200, 100, "Prs-coef.:", &cnt_pressure, -1.0e10, 1.0e10, 0.01, 5);
    wall_pack->end();
    wall_frame->end();

    // Ring パネルです。
    const int ring_panel_width = cnt_settings_column_width;
    const int ring_panel_height = 130;

    Fl_Group* ring_panel = create_framed_group_xy(0, 0, ring_panel_width - 50, ring_panel_height, "Ring");
    Fl_Pack* ring_pack = create_panel_pack_xy(10, 30, ring_panel_width - 50, ring_panel_height, true);

    create_float_input(200, 100, "Radius(A):", &cnt_ring_radius);
    create_float_input(200, 100, "F_load(eV/A):", &cnt_ring_fmax);
    create_float_input(200, 100, "Sharpness:", &cnt_ring_sharpness);
    ring_pack->end();
    ring_panel->end();

    cnt_settings_column2->end();
    cnt_settings_row1->end();

    // Which layer to show パネルです。
    const int which_layer_panel_width = cnt_setting_panel_width - 60;
    const int which_layer_panel_height = 70;

    FramedPanel* which_layer_panel = new FramedPanel(
        cnt_setting_panel_width, which_layer_panel_width, which_layer_panel_height, "Which layer to show", false);
    which_layer_panel->begin();

    create_check_button(60, "All", &show_cnt_all);
    create_check_button(60, "1st", &show_cnt_layer[0]);
    create_check_button(60, "2nd", &show_cnt_layer[1]);
    create_check_button(60, "3rd", &show_cnt_layer[2]);
    create_check_button(60, "4th", &show_cnt_layer[3]);
    create_check_button(60, "5th", &show_cnt_layer[4]);
    create_check_button(60, "6th", &show_cnt_layer[5]);
    which_layer_panel->begin();
    cnt_setting_panel->end();

    create_vertical_space(10);

    // Set parameters ウィンドウ 10 行目を作成します。
    Fl_Pack* setup_row10 = create_panel_pack(window_width, window_width - 100, 70, false);
    setup_row10->spacing(20);
    create_file_chooser_button(150, "Read SETDAT", pointer_cb, "SETDAT", CB_SETDAT_FB);
    create_process_button(150, "Write SETDAT", pointer_cb, WRITE_SETDAT_ID);
    create_window_close_button(150, "Close");
    setup_row10->end();

    create_vertical_space(10);

    end_window();
}

int SetParameterWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 何もしません。
            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            // 何もしません。
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}
