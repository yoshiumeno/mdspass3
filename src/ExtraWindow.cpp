#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include "ExtraWindow.h"
#include "FramedPanel.h"
#include "RolloutPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

ExtraWindow::ExtraWindow() : BaseWindow(550, 750, "Extra analysis")
{
    // 別ウィンドウのインスタンスを作成します。
    end_window();
    measure_window = std::make_unique<MeasureWindow>();
    begin_window();

    // Close ボタン以外の要素を含むパネルを作成します。
    Fl_Pack* panel1 = create_row(window_width, 700, 20);

    // 列のサイズを指定します。
    const int column_width = window_width / 2 - 10;
    const int column_heght = window_height;

    // パネル 1 の 1 列目です。
    Fl_Pack* column1 = create_column(column_width, column_heght);

    create_vertical_space(20);

    // Measure ボタンを作成します。
    create_window_open_button(80, "Measure", measure_window);

    // Force check パネルを作成します。
    const int force_check_width = column_width;
    const int force_check_height = 110;

    Fl_Group* force_check_frame = create_framed_group_xy(
        0, 0, force_check_width, force_check_height, "Force check");

    Fl_Pack* force_check_pack = create_panel_pack_xy(
        10, 30, force_check_width - 20, force_check_height, true);

    create_float_input(force_check_width - 60, 70, "disp:", &fcheck_disp, CB_FCHECK);
    create_vertical_space(10);
    create_process_button(80, "Check", control_cb, CB_FCHECK);

    force_check_pack->end();
    force_check_frame->end();

    // Phonon/DoS パネルを作成します。
    const int phonon_width = column_width;
    const int phonon_height = 260;

    Fl_Group* phonon_frame = create_framed_group_xy(
        0, 0, phonon_width, phonon_height, "Phonon/DoS");

    Fl_Pack* phonon_pack = create_panel_pack_xy(
        10, 30, phonon_width - 20, phonon_height, true);

    create_int_input(phonon_width - 40, 70, "# of replica:", &phonon_rep, CB_PHONON);
    create_dropdown_list(phonon_width - 40, 70, "k path", &phonon_kp, kpstring_list, MAXKP, CB_PHONON);
    create_int_input(phonon_width - 40, 70, "# of k-points:", &phonon_knum, CB_PHONON);
    create_vertical_space(10);
    create_process_button(80, "Calc (E-k)", pointer_cb, PHONON_CALC_ID);
    create_vertical_space(10);
    create_float_input(phonon_width - 40, 70, "Gauss width:", &dos_gauss_width);
    create_int_input(phonon_width - 40, 70, "DoS k-mesh", &dos_kmesh);
    create_vertical_space(10);
    create_process_button(80, "Calc (DoS)", pointer_cb, DOS_CALC_ID);

    phonon_pack->end();
    phonon_frame->end();

    // NEB パネルを作成します。
    const int neb_width = column_width;
    const int neb_height = 270;

    Fl_Group* neb_frame = create_framed_group_xy(
        0, 0, neb_width, neb_height, "NEB");
    Fl_Pack* neb_pack = create_panel_pack_xy(
        10, 30, neb_width - 20, neb_height - 30, true);

    create_check_button(neb_width - 40, "Init conf from file", &neb_init_read);
    create_int_input(neb_width - 40, 70, "# of nodes:", &neb_num, CB_NEB);
    create_int_input(neb_width - 40, 70, "Max iter:", &neb_ite, CB_NEB);
    create_float_input(neb_width - 40, 70, "tol 10^:", &neb_tol_fac, CB_NEB);
    create_check_button(neb_width - 40, "const spring", &neb_kconst, CB_NEB);
    create_float_input(neb_width - 40, 70, "k factor:", &neb_kfactor, CB_NEB);
    create_vertical_space(10);
    create_process_button(80, "Calc", pointer_cb, NEB_CALC_ID);
    create_vertical_space(10);
    create_int_spinner(neb_width - 40, 70, "Node to show", &neb_node, 0, 10000, CB_NEB_SHOW);

    neb_pack->end();
    neb_frame->end();

    column1->end();

    // パネル 1 の 2 列目です。
    Fl_Pack* column2 = create_column(column_width, column_heght);

    create_vertical_space(10);

    // CNT SHELL ANALYSIS パネルを作成します。
    const int cnt_shell_width = column_width;
    const int cnt_shell_height = 230;

    Fl_Group* cnt_shell_frame = create_framed_group_xy(
        0, 0, cnt_shell_width, cnt_shell_height, "CNT SHELL ANALYSIS");
    Fl_Pack* cnt_shell_pack = create_panel_pack_xy(
        10, 30, cnt_shell_width - 20, cnt_shell_height - 30, true);

    create_int_input(cnt_shell_width - 40, 70, "n:", &cntshell_n);
    create_float_input(cnt_shell_width - 40, 70, "d_mu:", &cntshell_dmu);
    create_float_input(cnt_shell_width - 40, 70, "d_nu:", &cntshell_dnu);
    create_vertical_space(10);
    create_process_button(80, "Dfm test", pointer_cb, CNTSHELL_TEST_ID);
    create_vertical_space(10);
    create_float_input(cnt_shell_width - 40, 70, "eps:", &cntshell_eps);
    create_vertical_space(10);
    create_process_button(80, "Calc", pointer_cb, CNTSHELL_CALC_ID);

    cnt_shell_pack->end();
    cnt_shell_frame->end();

    // Auto calc パネルを作成します。
    const int auto_calc_width = column_width;
    const int auto_calc_height = 135;

    Fl_Group* auto_calc_frame = create_framed_group_xy(
        0, 0, auto_calc_width, auto_calc_height, "Auto calc");
    Fl_Pack* auto_calc_pack = create_panel_pack_xy(
        10, 30, auto_calc_width - 20, auto_calc_height - 30, true);

    create_float_input(auto_calc_width - 40, 70, "val1:", &auto_val1);
    create_float_input(auto_calc_width - 40, 70, "val2:", &auto_val2);
    create_vertical_space(10);
    create_process_button(80, "Auto", pointer_cb, EXTRA_AUTO_ID);

    auto_calc_pack->end();
    auto_calc_frame->end();

    // Bookkeep table パネルを作成します。
    const int bookkeep_width = column_width;
    const int bookkeep_height = 100;

    Fl_Group* bookkeep_frame = create_framed_group_xy(
        0, 0, bookkeep_width, bookkeep_height, "Bookkeep table");
    Fl_Pack* bookkeep_pack = create_panel_pack_xy(
        10, 30, bookkeep_width - 20, bookkeep_height - 30, true);

    create_process_button(80, "Write", pointer_cb, EXTRA_BKWRITE_ID);
    create_process_button(80, "Read", pointer_cb, EXTRA_BKREAD_ID);

    bookkeep_pack->end();
    bookkeep_frame->end();

    column2->end();
    panel1->end();

    /*
    // Rigid_shift パネルを作成
    const int rigid_shift_width = column_width;
    //const int rigid_shift_width = 600; //instead of column_width
    const int rigid_shift_height = 440;
    RolloutPanel* rigid_shift_panel =
      //new RolloutPanel(window_width, rigid_shift_width, rigid_shift_height,
        new RolloutPanel(column_width, rigid_shift_width, rigid_shift_height,
                "Rigid shift (xy plane)", false);
    rigid_shift_panel->spacing(10);
    rigid_shift_panel->begin();

    // Rigit shift パネルの 1 列目です。
    //const int rigid_shift_column_width = rigid_shift_width / 3 - 20;
    const int rigid_shift_column_width = rigid_shift_width - 50;
    const int rigid_shift_column_height = rigid_shift_height;
    Fl_Pack* rigid_shift_column1 =
        create_column(rigid_shift_column_width, rigid_shift_column_height);

    create_vertical_space(10);

    create_float_input(rigid_shift_column_width, 80, "dx:", &rigid_shift_dx);
    create_float_input(rigid_shift_column_width, 80, "dy:", &rigid_shift_dy);
    create_float_input(rigid_shift_column_width, 80, "dz:", &rigid_shift_dz);

    ////rigid_shift_column1->end();

    // Rigit shift パネルの 2 列目です。
    ////Fl_Pack* rigid_shift_column2 =
    ////    create_column(rigid_shift_column_width, rigid_shift_column_height);
    create_vertical_space(10);
    create_float_input(rigid_shift_column_width, 60,
        "Separete at (z):", &rigid_shift_z);
    create_vertical_space(10);
    create_process_button(80, "Shift", control_cb, CB_RIGID_SHIFT);
    ////rigid_shift_column2->end();

    // Rigit shift パネルの 3 列目です。
    ////Fl_Pack* rigid_shift_column3 =
    ////    create_column(rigid_shift_column_width, rigid_shift_column_height);
    create_vertical_space(10);
    create_int_input(rigid_shift_column_width, 60,
        "times (x):", &rigid_shift_xtimes);
    create_int_input(rigid_shift_column_width, 60,
        "times (y):", &rigid_shift_ytimes);
    create_vertical_space(10);
    create_process_button(80, "Shift X-Y", control_cb, CB_RIGID_SHIFT_XY);
    ////rigid_shift_column3->end();
    rigid_shift_column1->end();

    rigid_shift_panel->end();

    column2->end();
    panel1->end();
    */

    // Close ボタンを作成します。
    create_window_close_button(80, "Close");

    end_window();
}

int ExtraWindow::handle(int event)
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
