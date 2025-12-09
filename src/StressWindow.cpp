#include "StressWindow.h"
#include "FramedPanel.h"
#include "RolloutPanel.h"
#include "myheader.h"

StressWindow::StressWindow() : BaseWindow(600, 600, "Stress")
{
    // テキストを追加します。
    new Fl_Box(0, 0, 200, 30, "Stress in MPa");

    // パネル 1 を作成します。
    const int panel1_width = 600;
    Fl_Pack* panel1 = create_panel_pack(window_width, panel1_width, 120, true);
    panel1->spacing(10);

    // パネル 1 の 1 行目です。
    Fl_Pack* panel1_row1 = create_row(panel1_width, 30, 10);
    create_horizontal_space(( window_width - panel1_width ) / 2);
    create_float_input(150, 100, "xx:", &strs_xx);
    create_float_input(150, 100, "yy:", &strs_yy);
    create_float_input(150, 100, "zz:", &strs_zz);
    panel1_row1->end();

    // パネル 1 の 2 行目です。
    Fl_Pack* panel1_row2 = create_row(panel1_width, 30, 10);
    create_horizontal_space(( window_width - panel1_width ) / 2);
    create_float_input(150, 100, "xy:", &strs_xy);
    create_float_input(150, 100, "yz:", &strs_yz);
    create_float_input(150, 100, "zx:", &strs_zx);
    panel1_row2->end();

    // パネル 1 の 3 行目です。
    Fl_Pack* panel1_row3 = create_row(panel1_width, 30, 10);
    create_horizontal_space(( window_width - panel1_width ) / 2);
    create_float_input(250, 100, "Stress check with de:", &eps_strschk);
    create_process_button(100, "Check", pointer_cb, STRSCHK_ID);
    panel1_row3->end();

    panel1->end();

    // PR Setting (Stress in MPa) パネルを作成します。
    const int pr_setting_width = 600;
    const int pr_setting_height = 220;
    RolloutPanel* pr_setting_panel = new RolloutPanel(
        window_width, pr_setting_width, pr_setting_height, "PR Setting (Stress in MPa)", true);
    pr_setting_panel->spacing(10);
    pr_setting_panel->begin();

    // PR Setting パネルの 1 行目です。
    create_float_input(200, 100, "PR mass x 10^:", &prmass_scale);

    // PR Setting パネルの 2 行目です。
    Fl_Pack* pr_setting_row2 = create_row(pr_setting_width, 30, 10);
    create_float_input(150, 100, "xx:", &strs_set_xx, STRSSET_ID, pointer_cb);
    create_float_input(150, 100, "yy:", &strs_set_yy, STRSSET_ID, pointer_cb);
    create_float_input(150, 100, "zz:", &strs_set_zz, STRSSET_ID, pointer_cb);
    pr_setting_row2->end();

    // PR Setting パネルの 3 行目です。
    Fl_Pack* pr_setting_row3 = create_row(pr_setting_width, 30, 10);
    create_float_input(150, 100, "xy:", &strs_set_xy, STRSSET_ID, pointer_cb);
    create_float_input(150, 100, "yz:", &strs_set_yz, STRSSET_ID, pointer_cb);
    create_float_input(150, 100, "zx:", &strs_set_zx, STRSSET_ID, pointer_cb);
    pr_setting_row3->end();

    // PR Setting パネルの 4 行目です。
    Fl_Pack* pr_setting_row4 = create_row(pr_setting_width, 30, 10);
    create_float_input(180, 100, "Damper:", &prdamper_val1);
    create_float_input(180, 100, "Hv Limit", &strs_set_zx);
    pr_setting_row4->end();

    // PR Setting パネルを終了します。
    pr_setting_panel->end();

    // Cell constraint パネルを作成します。
    const int cell_constraint_width = 500;
    const int cell_constraint_height = 80;
    RolloutPanel* cell_constraint_panel = new RolloutPanel(
        window_width, cell_constraint_width, cell_constraint_height, "Cell constraint", false);
    cell_constraint_panel->spacing(20);
    cell_constraint_panel->begin();
    create_check_button(50, "xx", &cellfix_xx, CB_CELLFIX_ID);
    create_check_button(50, "yy", &cellfix_yy, CB_CELLFIX_ID);
    create_check_button(50, "zz", &cellfix_zz, CB_CELLFIX_ID);
    create_check_button(50, "xy", &cellfix_xy, CB_CELLFIX_ID);
    create_check_button(50, "yz", &cellfix_yz, CB_CELLFIX_ID);
    create_check_button(50, "zx", &cellfix_zx, CB_CELLFIX_ID);
    cell_constraint_panel->end();

    // CNT corrugation パネルを作成します。
    const int cnt_corrugation_width = 600;
    const int cnt_corrugation_height = 120;
    RolloutPanel* cnt_corrugation_panel = new RolloutPanel(
        window_width, cnt_corrugation_width, cnt_corrugation_height, "CNT corrugation", false);
    cnt_corrugation_panel->spacing(10);
    cnt_corrugation_panel->begin();

    // CNT corrugation パネル の 1 列目です。
    const int cnt_corrugation_column1_width = 150;
    Fl_Pack* cnt_corrugation_column1 = create_column(
        cnt_corrugation_column1_width, cnt_corrugation_height);
    create_vertical_space(10);
    create_float_input(150, 80, "Load(nN):", &cnt_pressure_ftot);
    cnt_corrugation_column1->end();

    // CNT corrugation パネル の 2 列目です。
    const int cnt_corrugation_column2_width = 260;
    Fl_Pack* cnt_corrugation_column2 = create_column(
        cnt_corrugation_column2_width, cnt_corrugation_height);
    create_vertical_space(10);
    create_float_input(260, 80, "Pressure (GPa):", &cnt_pressure_gpa);
    create_float_input(260, 80, "Outermost CNT area (A^2):", &cylinder_side_area_f);
    cnt_corrugation_column2->end();

    // CNT corrugation パネル の 3 列目です。
    const int cnt_corrugation_column3_width = 130;
    Fl_Pack* cnt_corrugation_column3 = create_column(
        cnt_corrugation_column3_width, cnt_corrugation_height);
    create_vertical_space(10);
    create_float_input(110, 80, "or:", &cnt_pressure_gpa2);
    create_float_input(110, 80, "or:", &cylinder_side_area0_f);
    cnt_corrugation_column3->end();

    cnt_corrugation_panel->end();

    // Auxiliary constraint パネルを作成します。
    const int auxiliary_constraint_width = 300;
    const int auxiliary_constraint_height = 80;
    RolloutPanel* auxiliary_constraint_panel = new RolloutPanel(
        window_width, auxiliary_constraint_width, auxiliary_constraint_height, "Auxiliary constraint", true);
    auxiliary_constraint_panel->begin();
    create_float_input(250, 80, "yz pln load(nN):", &yzplane_punch_ftot);
    auxiliary_constraint_panel->end();

    create_vertical_space(10);
    Fl_Pack* button_column = create_column(200, 80, 10);
    create_process_button(180, "Calc elastic coeff.", control_cb, CB_ELASC);
    create_process_button(180, "Close", close_window_callback);
    create_vertical_space(10);
    button_column->end();

    end_window();
}

int StressWindow::handle(int event)
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
