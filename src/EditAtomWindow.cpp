#include "EditAtomWindow.h"
#include "FramedPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

void set_atom_color();

EditAtomWindow::EditAtomWindow()
    : BaseWindow(700, 660, "Edit atom")
{
    // Edit atom ウィンドウの 1 行目です。
    // パネル 1 を作成します。
    const int panel1_width = 600;
    const int panel1_height = 370;
    FramedPanel* panel1 = new FramedPanel(
        window_width, panel1_width, panel1_height, "", true);
    panel1->begin();

    // パネル 1 の 1 行目です。
    const int panel1_row1_width = panel1_width;
    const int panel1_row1_height = 180;
    Fl_Pack* panel1_row1 = create_row(panel1_row1_width, panel1_row1_height, 30);

    // パネル 1 の 1 行目の 1 列目です。
    const int panel1_row1_column_width = panel1_row1_width / 3 - 50;
    const int panel1_row1_column_height = panel1_row1_height;
    Fl_Pack* panel1_row1_column1 = create_column(panel1_row1_column_width, panel1_row1_height);
    create_vertical_space(10);

    //create_int_spinner(panel1_row1_column_width, 80, "atom #:", &iatom_pick, CB_EDITATOM);
    create_int_spinner(panel1_row1_column_width, 80, "atom #:", &iatom_pick, 1, atom.natom, CB_EDITATOM); //YU2024.08
    create_text_input(panel1_row1_column_width, 80, "Atom:", config_atom);
    create_float_input(panel1_row1_column_width, 80, "x (A):", &atomrx);
    create_float_input(panel1_row1_column_width, 80, "y (A):", &atomry);
    create_float_input(panel1_row1_column_width, 80, "z (A):", &atomrz);

    panel1_row1_column1->end();

    // パネル 1 の 1 行目の 2 列目です。
    Fl_Pack* panel1_row1_column2 = create_column(panel1_row1_column_width, panel1_row1_height);
    create_vertical_space(10);
    create_process_button(100, "Remove atom", pointer_cb, REMOVEATOM_DO_ID);
    create_process_button(100, "Add atom", pointer_cb, ADDATOM_DO_ID);
    create_process_button(100, "Move/Change", pointer_cb, MOVEATOM_DO_ID);
    create_process_button(100, "Shift one", pointer_cb, SHIFTONE_DO_ID);
    create_process_button(100, "Shift all", pointer_cb, SHIFTALL_DO_ID);
    panel1_row1_column2->end();

    // パネル 1 の 1 行目の 3 列目です。
    Fl_Pack* panel1_row1_column3 = create_column(panel1_row1_column_width, panel1_row1_height);

    create_vertical_space(10);
    create_float_input(panel1_row1_column_width, 100, "Rot x:", &rotatomx);
    create_float_input(panel1_row1_column_width, 100, "Rot y:", &rotatomy);
    create_float_input(panel1_row1_column_width, 100, "Rot z:", &rotatomz);
    create_process_button(100, "Rotate all", pointer_cb, ROTATEALL_DO_ID);

    panel1_row1_column3->end();
    panel1_row1->end();

    // Partial shift/rotation パネルを作成します。
    const int partial_shift_width = 400;
    const int partial_shift_height = 140;
    ::Fl_Group* partial_shift_frame = create_framed_group_xy(
        0, 0, partial_shift_width, partial_shift_height, "Partial shift/rotation");

    Fl_Pack* partial_shift_pack = create_panel_pack_xy(
        0, 30, partial_shift_width, partial_shift_height - 30, false);

    // Partial shift/rotation パネル の 1 列目です。
    const int partial_shift_column_width = 150;
    create_horizontal_space(80);

    Fl_Pack* partial_shift_column1 = create_column(partial_shift_column_width, partial_shift_height);
    create_vertical_space(10);
    create_int_spinner(partial_shift_column_width, 100, "From:", &iatom_from, -10000, 10000, CB_EDITATOM);
    create_int_spinner(partial_shift_column_width, 100, "To:", &iatom_to, -10000, 10000, CB_EDITATOM);
    partial_shift_column1->end();

    create_horizontal_space(20);

    // Partial shift/rotation パネル の 2 列目です。
    Fl_Pack* partial_shift_column2 = create_column(partial_shift_column_width, partial_shift_height);
    create_vertical_space(10);
    create_process_button(80, "Shift part", pointer_cb, SHIFTPART_DO_ID);
    create_process_button(80, "Rotate part", pointer_cb, ROTATEPART_DO_ID);
    create_process_button(80, "Remove part", pointer_cb, REMOVEPART_DO_ID);
    partial_shift_column2->end();

    partial_shift_pack->end();
    partial_shift_frame->end();

    panel1->end();

    /*
    // Edit atom ウィンドウ中程に追加（-->やっぱりやめた）
    // パネル２
    const int panel2_width = 600;
    const int panel2_height = 250;
    FramedPanel* panel2 = new FramedPanel(
        window_width, panel2_width, panel2_height, "Atom grouping", true);
    panel2->begin();
    Fl_Pack* panel2_row1 = create_row(600,250);
    create_horizontal_space(50);
    Fl_Pack* panel2_column1 = create_column(150,250);
    create_int_spinner(150, 100, "Group ID", &group_num, 0, atom.natom);
    create_radio_button(150, "Range in x", &group_range_axis, 0);
    create_radio_button(150, "Range in y", &group_range_axis, 1);
    create_radio_button(150, "Range in z", &group_range_axis, 2);
    create_float_input(150, 80, "max", &group_range_max);
    create_float_input(150, 80, "min", &group_range_min);
    create_process_button(80, "Set group", pointer_cb, SETGROUP_DO_ID);
    create_process_button(80, "Reset group", pointer_cb, RESETGROUP_DO_ID);
    panel2_column1->end();
    create_horizontal_space(50);
    Fl_Pack* panel2_column2 = create_column(150,250);
    create_check_button(50, "Visualize group", &if_group_num_show, CB_SHOW_GROUP);
    create_int_spinner(150, 80, "#:", &group_num_show, 0, MAXGROUP, CB_SHOW_GROUP);
    panel2_column2->end();
    panel2_row1->end();
    panel2->end();
    */

    // Edit atom ウィンドウの 2 行目です。
    // Polyhedron パネルを作成します。
    const int polyhedron_width = 650;
    const int polyhedron_height = 180;
    FramedPanel* polyhedron_panel = new FramedPanel(
        window_width, polyhedron_width, polyhedron_height, "Polyhedron", false);
    polyhedron_panel->begin();

    // Polyhedron パネルの 1 列目です。
    const int polyhedron_column_width = polyhedron_width / 4 - 20;
    Fl_Pack* polyhedron_column1 = create_column(polyhedron_column_width, polyhedron_height);
    create_vertical_space(10);
    create_radio_button(100, "Tetrahedron", &polyhedtype, 0);
    create_radio_button(100, "Hexahedron", &polyhedtype, 1);
    create_radio_button(100, "Octahedron", &polyhedtype, 2);
    create_radio_button(100, "Dodecahedron", &polyhedtype, 3);
    create_radio_button(100, "Icosahedron", &polyhedtype, 4);
    polyhedron_column1->end();

    create_horizontal_space(10);

    // Polyhedron パネルの 2 列目です。
    Fl_Pack* polyhedron_column2 = create_column(polyhedron_column_width, polyhedron_height);
    create_vertical_space(10);
    create_text_input(polyhedron_column_width, 80, "Atom:", config_atom);
    polyhedron_column2->end();

    create_horizontal_space(10);

    // Polyhedron パネルの 3 列目です。
    Fl_Pack* polyhedron_column3 = create_column(polyhedron_column_width, polyhedron_height);
    create_vertical_space(10);
    create_float_input(polyhedron_column_width, 80, "Size:", &polyhedsize);
    polyhedron_column3->end();

    create_horizontal_space(10);

    // Polyhedron パネルの 4 列目です。
    Fl_Pack* polyhedron_column4 = create_column(polyhedron_column_width, polyhedron_height);
    create_vertical_space(10);
    create_process_button(100, "Add Polyhedron", pointer_cb, ADDPOLYHED_DO_ID);
    polyhedron_column4->end();

    polyhedron_panel->end();

    // ボタンを配置します。
    Fl_Pack* button_row = create_row(window_width, 25, 20);
    create_horizontal_space(( window_width - 250 ) / 2);
    create_process_button(100, "Pack all in cell", pointer_cb, PACKALL_DO_ID);
    create_window_close_button(100, "Close");
    button_row->end();
}

int EditAtomWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            editatom_mode = 1;
            return 1;

            // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            editatom_mode = 0;
            set_atom_color();
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}
