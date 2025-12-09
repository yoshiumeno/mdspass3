#include "MainWindow.h"
#include "FramedPanel.h"
#include "RolloutPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

MainWindow::MainWindow() : BaseWindow(600, 600, "Control")
{
    end_window();
    // メイン ウィンドウから呼び出す他のウィンドウのインスタンスを作成します。
    edit_element_window = std::make_unique<EditElementWindow>();
    edit_element_xz_plane_window = std::make_unique<EditElementXZPlaneWindow>();
    set_parameter_window = std::make_unique<SetParameterWindow>();
    stress_window = std::make_unique<StressWindow>();
    extra_window = std::make_unique<ExtraWindow>();
    create_config_window = std::make_unique<CreateConfigWindow>();
    potential_arg_window = std::make_unique<PotentialArgWindow>();
    md_viewer_window = std::make_unique<MdViewerWindow>();
    open_file_window = std::make_unique<OpenFileWindow>();
    quit_window = std::make_unique<QuitWindow>();

    // MD Viewer ウィンドウを表示します。
    md_viewer_window->show();

    // メイン ウィンドウの作成を開始します。
    begin_window();

    // テキストを表示します。
    Fl_Box* title_text = new Fl_Box(0, 0, 280, 25, "MDSPASS ver.3.0");

    // 区切り線です。
    Fl_Box* separator = new Fl_Box(0, 0, 480, 2);
    separator->box(FL_FLAT_BOX);
    separator->color(FL_BLACK);

    // 1 行目です。
    // PBC に関するチェック ボックスを作成します。
    Fl_Pack* row1 = create_row(window_width, 30, 10);
    create_horizontal_space(( window_width - 300 ) / 2 - 20);
    create_check_button(90, "PBC x", &cell.pbcx);
    create_check_button(90, "PBC y", &cell.pbcy);
    create_check_button(90, "PBC z", &cell.pbcz);
    row1->end();

    create_vertical_space(10);

    // 2 行目です。
    Fl_Pack* row2 = create_row(window_width, 30, 20);
    Fl_Choice* potential_choice = create_dropdown_list_for_potential_set(180, 100, "Potential:", &ipottype);
    for (int i = 0; i < MAXPOTTYPE; ++i)
    {
        potential_choice->add(atom.potstring_list[i]);
    }
    potential_choice->value(ipottype);

    create_text_input(150, 100, "ARG:", atom.potential_arg);
    create_window_open_button(130, "ARG/FILE", potential_arg_window);
    row2->end();

    create_vertical_space(10);

    // 3 行目です。パネルを表示します。
    Fl_Pack* row3 = create_panel_pack(window_width, 520, 80, false);
    row3->spacing(10);

    create_horizontal_space(50);

    Fl_Pack* row3_column1 = create_column(140, 80);
    create_window_open_button(100, "Set param", set_parameter_window);
    create_process_button(100, "Calc", pointer_cb, CALC_ID);
    row3_column1->end();

    Fl_Pack* row3_column2 = create_column(140, 80);
    create_process_button(100, "MD on/off", pointer_cb, MDSWITCH_ID);
    create_window_open_button(100, "Stress", stress_window);
    row3_column2->end();

    Fl_Pack* row3_column3 = create_column(140, 80);
    create_float_input(140, 60, "Temp set:", &temp_set);
    create_window_open_button(100, "Extra", extra_window);
    row3_column3->end();
    row3->end();

    // 4 行目です。QC settings パネルを作成します。
    const int qc_setting_width = 500;
    const int qc_setting_height = 150;
    RolloutPanel* qc_setting_panel = new RolloutPanel(
        window_width, qc_setting_width, qc_setting_height, "QC settings", true);
    qc_setting_panel->begin();

    Fl_Pack* qc_setting_row1 = create_row(qc_setting_width, 30, 10);
    create_check_button(200, "QC", &atom.QC);
    create_check_button(200, "Show only QC element", &show_only_elem);
    qc_setting_row1->end();

    Fl_Pack* qc_setting_row2 = create_row(qc_setting_width, 30, 10);
    create_window_open_button(200, "Edit elem", edit_element_window);
    create_window_open_button(200, "Edit elem xz", edit_element_xz_plane_window);
    qc_setting_row2->end();

    Fl_Pack* qc_setting_row3 = create_row(qc_setting_width, 30, 10);
    create_file_chooser_button(200, "Read qcelement", pointer_cb, "", CB_QCELM_FB);
    create_process_button(200, "Write qcelement", pointer_cb, WRITEQCELEMENT_ID);
    qc_setting_row3->end();

    qc_setting_panel->end();

    // 5 行目です。
    Fl_Pack* row5 = create_row(window_width, 50, 10);

    create_horizontal_space(100);

    Fl_Pack* row5_column1 = create_column(300, 80);
    create_dropdown_list(300, 220, "Algorithm:", &ensemble, ensemble_list, MAXENSTYPE, CB_ENSEMBLE);
    row5_column1->end();

    Fl_Pack* row5_column2 = create_column(120, 80);
    create_check_button(120, "No translation", &notrans);
    create_check_button(120, "Keep in cell", &incell);
    row5_column2->end();

    row5->end();

    // 6 行目です。Instability analysis パネルを作成します。
    const int instability_analysis_width = 500;
    const int instability_analysis_height = 230;
    RolloutPanel* instability_analysis_panel = new RolloutPanel(
        window_width, instability_analysis_width, instability_analysis_height, "Instability analysis", false);
    instability_analysis_panel->begin();

    Fl_Pack* instability_analysis_column1 = create_column(250, instability_analysis_height);
    create_vertical_space(10);
    create_check_button(120, "Evector", &ievec);
    create_int_spinner(200, 60, "Center atom (inst):", &atom.instcenter, 0, 10000, 2);
    create_radio_button(200, "instability_atomcell", &inst_mode, 0);
    create_radio_button(200, "instability_atom", &inst_mode, 1);
    create_radio_button(200, "instability_atom_noremove", &inst_mode, 2);
    create_check_button(120, "Read Hessian", &hessian_read);
    create_check_button(120, "Write Hessian", &hessian_write);
    instability_analysis_column1->end();

    Fl_Pack* instability_analysis_column2 = create_column(200, instability_analysis_height);
    create_vertical_space(10);
    create_int_spinner(200, 100, "Mode:", &ievec_num, 0, 10000, CB_INST);
    create_float_spinner(200, 100, "Length:", &evec_len, -10000, 10000, 0.5, CB_INST);
    create_float_spinner(200, 100, "Eig val:", &eigval, -10000, 10000, 0.5, CB_INST);
    create_vertical_space(10);
    create_process_button(150, "Inst analysis", pointer_cb, INST_ID);
    instability_analysis_column2->end();

    instability_analysis_panel->end();

    create_vertical_space(10);

    // 7 行目です。Capture 関連の GUI を作成します。
    Fl_Pack* row7 = create_row(window_width, 30, 10);
    create_horizontal_space(100);
    create_button_with_callback(100, "Capture", capture_callback);
    create_int_input(130, 60, "Cap-file#", &capture_count);
    create_process_button(100, "Reset view", control_cb, CB_RESET_VIEW);
    row7->end();

    create_vertical_space(10);

    // 8 行目です。Status パネルを作成します。
    const int status_width = 550;
    const int status_height = 170;
    RolloutPanel* status_panel = new RolloutPanel(
        window_width, status_width, status_height, "Status", false);
    // Status パネルを表示状態にします。
    status_panel->toggle();
    status_panel->spacing(10);
    status_panel->begin();

    Fl_Pack* status_column1 = create_column(280, status_height);
    create_vertical_space(10);
    create_int_input(250, 140, "Step:", &istep);
    create_float_input(250, 140, "E_P(eV/atm):", &epotatom);
    create_float_input(250, 140, "Temp(K):", &tempc);
    create_float_input(250, 140, "Fmax(eV/A):", &f_max);
    status_column1->end();

    Fl_Pack* status_column2 = create_column(200, status_height);
    create_vertical_space(10);
    create_float_input(200, 140, "Cell(A):", &cellx, STATUS_LX, pointer_cb);
    create_float_input(200, 140, "", &celly, STATUS_LY, pointer_cb);
    create_float_input(200, 140, "", &cellz, STATUS_LZ, pointer_cb);
    create_float_input(200, 140, "dt (fs):", &dtm, STATUS_DT, pointer_cb);
    status_column2->end();

    status_panel->end();

    // 9 行目です。File control and config createtion パネルを作成します。
    const int file_width = 500;
    const int file_height = 100;
    RolloutPanel* file_panel = new RolloutPanel(
        window_width, file_width, file_height, "File control and config creation", true);
    // File control and config creation パネルを表示状態にします。
    file_panel->toggle();
    file_panel->spacing(10);
    file_panel->begin();

    Fl_Pack* file_row1 = create_row(file_width, 30, 10);

    create_window_open_button(150, "Read Config", open_file_window);
    create_process_button(150, "Write Config", pointer_cb, WRITECONFIG_ID);

    create_window_open_button(150, "Create Config", create_config_window);

    file_row1->end();

    file_panel->end();

    create_vertical_space(10);

    // 10 行目です。Reset、Close ボタンを配置するパネルを作成します。
    Fl_Pack* row10 = create_row(window_width, 30, 10);
    create_horizontal_space(( window_width - 220 ) / 2);
    create_process_button(100, "Reset", pointer_cb, RESET_ID);
    //create_window_close_button(100, "Close");
    create_window_open_button(100, "Quit", quit_window);

    row10->end();

    end_window();
}

int MainWindow::handle(int event)
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
