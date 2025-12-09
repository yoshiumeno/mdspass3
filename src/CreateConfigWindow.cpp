#include "CreateConfigWindow.h"
#include "FramedPanel.h"
#include "RolloutPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

// 関数のプロトタイプ宣言です。
void set_atom_color();

CreateConfigWindow::CreateConfigWindow()
    : BaseWindow(600, 850, "Create Config")
{
    end_window();
    // Edit atom ウィンドウのインスタンスを作成します。
    edit_atom_window = std::make_unique<EditAtomWindow>();
    begin_window();

    // Edit atom ウィンドウの 1 行目のパネルを作成します。
    const int panel1_width = 580;
    const int panel1_height = 380;
    FramedPanel* panel1 =
        new FramedPanel(window_width, panel1_width, panel1_height, "", false);
    panel1->begin();

    // パネル 1 の 1 列目です。
    const int panel1_column_width = panel1_width / 2 - 10;
    const int panel1_column_height = panel1_height;
    Fl_Pack* panel1_column1 =
        create_column(panel1_column_width, panel1_column_height);

    create_vertical_space(10);

    create_text_input(panel1_column_width - 30, 100, "Atom:", config_atom);
    create_text_input(panel1_column_width - 30, 100, "Atom(2):", config_atom2);
    // ドロップダウン リストを作成します。
    Fl_Choice* potential_choice = create_dropdown_list_for_potential_set(
        panel1_column_width - 30, 100, "Potential:", &ipottype);
    for (int i = 0; i < MAXPOTTYPE; ++i)
    {
        potential_choice->add(atom.potstring_list[i]);
    }
    potential_choice->value(ipottype);
    create_text_input(panel1_column_width - 30, 100, "ARG:", atom.potential_arg);

    // ラジオ ボタンを作成します。
    create_radio_button(panel1_column_width - 30, "FCC", &config_type, 0, 0);
    create_radio_button(panel1_column_width - 30, "BCC", &config_type, 1, 0);
    create_radio_button(panel1_column_width - 30, "Diamond", &config_type, 2, 0);
    create_radio_button(panel1_column_width - 30, "Wurtzite", &config_type, 3, 0);
    create_radio_button(panel1_column_width - 30, "Nanotube", &config_type, 4, 0);
    create_radio_button(panel1_column_width - 30, "2-D Triangle", &config_type, 5,
        0);
    create_radio_button(panel1_column_width - 30, "Graphene", &config_type, 6, 0);

    // スピナーを作成します。
    create_float_spinner(panel1_column_width - 30, 80,
        "NT Rotation (z) [deg]:", &rotz, 0.0, 360.0, 0.5, 0);
    create_float_spinner(panel1_column_width - 30, 80,
        "NT Shift (z) [%]:", &shiftz, 0.0, 100.0, 0.5, 0);
    panel1_column1->end();

    create_vertical_space(10);

    // パネル 1 の 2 列目です。
    Fl_Pack* panel1_column2 =
        create_column(panel1_column_width, panel1_column_height);

    create_vertical_space(10);

    // スピナーを作成します。
    create_int_spinner(panel1_column_width - 30, 60, "# of rep in x:", &irepx, 1,
        1000, 0);
    create_int_spinner(panel1_column_width - 30, 60, "# of rep in y:", &irepy, 1,
        1000, 0);
    create_int_spinner(panel1_column_width - 30, 60, "# of rep in z:", &irepz, 1,
        1000, 0);
    create_int_spinner(panel1_column_width - 30, 60, "m of (m, n):", &icntm, 1,
        1000, 0);
    create_int_spinner(panel1_column_width - 30, 60, "n of (m, n):", &icntn, 1,
        1000, 0);
    create_float_spinner(panel1_column_width - 30, 60,
        "NT cell size (x/y) [ang]:", &cscnt, 10.0f, 200.0, 0.5,
        0);
    create_float_spinner(panel1_column_width - 30, 60,
        "Lattice const a [ang]:", &alat, 0.0f, 100.0, 0.5, 0);
    create_float_spinner(panel1_column_width - 30, 60,
        "Lattice const c [ang]:", &clat, 0.0f, 100.0, 0.5, 0);

    panel1_column2->end();

    panel1->end();

    // Create lattice パネルを作成します。
    const int create_lattice_width = 600;
    const int create_lattice_height = 160;
    RolloutPanel* create_lattice_panel =
        new RolloutPanel(window_width, create_lattice_width,
            create_lattice_height, "Create lattice", false);
    create_lattice_panel->spacing(20);
    create_lattice_panel->begin();

    // Creaqet lattice パネルの 1 列目です。
    const int create_lattice_column1_width = 150;
    const int create_lattice_column1_height = create_lattice_height;
    Fl_Pack* create_lattice_column1 = create_column(
        create_lattice_column1_width, create_lattice_column1_height);
    create_vertical_space(10);
    create_radio_button(create_lattice_column1_width - 10, "Rhombohedral",
        &lattice_type, 0, 0);
    create_lattice_column1->end();

    // Create lattice パネルの 2 列目です。
    const int create_lattice_column2_width = 350;
    const int create_lattice_column2_height = create_lattice_height;
    Fl_Pack* create_lattice_column2 = create_column(
        create_lattice_column2_width, create_lattice_column2_height);
    create_vertical_space(10);
    create_float_input(create_lattice_column2_width - 30, 100,
        "Lattice const. (a) [ang]:", &alat);
    create_float_input(create_lattice_column2_width - 30, 100,
        "Lattice const. (c) [ang]:", &clat);
    create_float_input(create_lattice_column2_width - 30, 100,
        "Angle 1 [deg]:", &alplat1);
    create_vertical_space(10);
    create_process_button(80, "Create Lattice", control_cb, CB_LATTICE);

    create_lattice_column2->end();

    create_lattice_panel->end();

    // Miller index (cubic cell only) パネルを作成します。
    const int miller_index_width = 400;
    const int miller_index_height = 120;
    RolloutPanel* miller_index_panel =
        new RolloutPanel(window_width, miller_index_width, miller_index_height,
            "Miller index (cubic cell only)", true);
    miller_index_panel->begin();

    // Miller index パネルの 1 行目です。
    const int miller_index_row_width = miller_index_width;
    const int miller_index_row_height = 25;
    Fl_Pack* miller_index_row1 =
        create_row(miller_index_row_width, miller_index_row_height, 20);

    create_int_input(80, 60, "1:", &milX1);
    create_int_input(80, 60, "", &milX2);
    create_int_input(80, 60, "", &milX3);

    miller_index_row1->end();

    // Miller index パネルの 2 行目です。
    Fl_Pack* miller_index_row2 =
        create_row(miller_index_row_width, miller_index_row_height, 20);
    create_int_input(80, 60, "2:", &milY1);
    create_int_input(80, 60, "", &milY2);
    create_int_input(80, 60, "", &milY3);
    miller_index_row2->end();

    // Miller index パネルの 3 行目です。
    Fl_Pack* miller_index_row3 =
        create_row(miller_index_row_width, miller_index_row_height, 20);
    create_int_input(80, 60, "3:", &milZ1);
    create_int_input(80, 60, "", &milZ2);
    create_int_input(80, 60, "", &milZ3);
    miller_index_row3->end();

    miller_index_panel->end();

    // Rotation パネルを作成します。
    const int rotation_width = 500;
    const int rotation_height = 120;
    RolloutPanel* rotation_panel = new RolloutPanel(
        window_width, rotation_width, rotation_height, "Rotation", true);
    rotation_panel->spacing(10);
    rotation_panel->begin();

    // Rotation パネルの 1 行目です。
    Fl_Pack* rotation_row1 = create_row(rotation_width, 25, 20);

    create_float_input(120, 100, "X:", &cellrot_x);
    create_float_input(120, 100, "Y:", &cellrot_y);
    create_float_input(120, 100, "Z:", &cellrot_z);

    rotation_row1->end();

    // Rotation パネルの 2 行目です。
    Fl_Pack* rotation_row2 = create_row(rotation_width, 25, 30);

    create_process_button(80, "Rotate", control_cb, CB_ROTATE);
    create_process_button(80, "Reverse", control_cb, CB_ROTATE_REV);

    rotation_row2->end();
    rotation_panel->end();

    // Cell dimension パネルを作成します。
    const int cell_dimension_width = 400;
    const int cell_dimension_height = 180;
    RolloutPanel* cell_dimension_panel =
        new RolloutPanel(window_width, cell_dimension_width,
            cell_dimension_height, "Cell dimension", true);
    cell_dimension_panel->begin();

    // Cell dimension パネル の 1 行目です。
    Fl_Pack* cell_dimension_row1 = create_row(cell_dimension_width, 30);

    create_float_input(80, 60, "1 ", &cell1x, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell1y, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell1z, CB_CELL_DIM);

    cell_dimension_row1->end();

    // Cell dimension パネル の 2 行目です。
    Fl_Pack* cell_dimension_row2 = create_row(cell_dimension_width, 30);

    create_float_input(80, 60, "2 ", &cell2x, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell2y, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell2z, CB_CELL_DIM);

    cell_dimension_row2->end();

    // Cell dimension パネル の 3 行目です。
    Fl_Pack* cell_dimension_row3 = create_row(cell_dimension_width, 30);

    create_float_input(80, 60, "3 ", &cell3x, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell3y, CB_CELL_DIM);
    create_float_input(80, 60, "", &cell3z, CB_CELL_DIM);

    cell_dimension_row3->end();

    create_vertical_space(10);

    // Cell dimension パネル の 4 行目です。
    Fl_Pack* cell_dimension_row4 = create_row(cell_dimension_width, 30);
    create_check_button(80, "Flx atoms", &ifix_atoms);
    cell_dimension_row4->end();
    cell_dimension_panel->end();

    // Atom grouping パネル作成
    const int atomgroup_width = 550;
    const int atomgroup_height = 250;
    RolloutPanel* atomgroup_panel =new RolloutPanel(window_width, atomgroup_width,
        atomgroup_height, "Atom Grouping", true);
    atomgroup_panel->begin();
    Fl_Pack* atomgroup_row1 = create_row(600,250);
    create_horizontal_space(50);
    Fl_Pack* atomgroup_column1 = create_column(150,250);
    create_int_spinner(150, 100, "Group ID", &group_num, 0, atom.natom);
    create_radio_button(150, "Range in x", &group_range_axis, 0);
    create_radio_button(150, "Range in y", &group_range_axis, 1);
    create_radio_button(150, "Range in z", &group_range_axis, 2);
    create_float_input(150, 80, "max", &group_range_max);
    create_float_input(150, 80, "min", &group_range_min);
    create_process_button(80, "Set group", pointer_cb, SETGROUP_DO_ID);
    create_process_button(80, "Reset group", pointer_cb, RESETGROUP_DO_ID);
    atomgroup_column1->end();
    create_horizontal_space(50);
    Fl_Pack* atomgroup_column2 = create_column(150,250);
    create_check_button(50, "Visualize group", &if_group_num_show, CB_SHOW_GROUP);
    create_int_spinner(150, 80, "#:", &group_num_show, 0, MAXGROUP, CB_SHOW_GROUP);
    atomgroup_column2->end();
    atomgroup_row1->end();
    atomgroup_panel->end();

    // Slicer パネルを作成します。
    const int slicer_width = 550;
    const int slicer_height = 440;
    RolloutPanel* slicer_panel = new RolloutPanel(window_width, slicer_width,
        slicer_height, "Slice/Trim atoms", true); //SlicerをSlice/Trimに変更
    slicer_panel->begin();

    // YU2025.01
    Fl_Group* slicer_frame = create_framed_group_xy(
        0, 0, slicer_width, 140, "Slicer");
    Fl_Pack* slicer_pack = create_panel_pack_xy(
        10, 30, slicer_width, 140, true);

    // Slicer パネルの 1 行目です。
    Fl_Pack* slicer_row1 = create_row(slicer_width, 25, 30);
    create_float_input(130, 80, "1 min:", &slice_1min);
    create_float_input(130, 80, "2 min:", &slice_2min);
    create_float_input(130, 80, "3 min:", &slice_3min);
    slicer_row1->end();

    // Slicer パネルの 2 行目です。
    Fl_Pack* slicer_row2 = create_row(slicer_width, 25, 30);
    create_float_input(130, 80, "1 max:", &slice_1max);
    create_float_input(130, 80, "2 max:", &slice_2max);
    create_float_input(130, 80, "3 max:", &slice_3max);
    slicer_row2->end();

    create_vertical_space(10);

    // Slicer パネルの 3 行目です。
    Fl_Pack* slicer_row3 = create_row(slicer_width, 30);
    create_process_button(180, "Slice", control_cb, CB_SLICE);
    slicer_row3->end();

    slicer_pack->end();
    slicer_frame->end();


    Fl_Group* trimmer_frame = create_framed_group_xy(0, 0, slicer_width, 160, "Trimmer");
    Fl_Pack* trimmer_pack = create_panel_pack_xy(10, 30, slicer_width, 160, true);
    // Trim パネル
    Fl_Pack* trimmer_row1 = create_row(slicer_width, 80, 10);
    Fl_Pack* trimmer_column1 = create_column(slicer_width/3,80);
    create_radio_button(50, "Cylinder", &trim_mode, 0);
    create_radio_button(50, "C-pillar", &trim_mode, 1);
    create_radio_button(50, "Hole", &trim_mode, 2);
    trimmer_column1->end();
    Fl_Pack* trimmer_column2 = create_column(slicer_width*2/3,80);
    Fl_Pack* trimmer_row11 = create_row(slicer_width, 25, 10);
    create_int_spinner(115, 80, "axis:", &trim_cylinder_axis, 1, 3);
    create_float_input(170, 80, "diameter(A):", &trim_cylinder_diameter);
    trimmer_row11->end();
    create_float_input(190, 80, "diameter 2(A):", &trim_cylinder_diameter2);
    trimmer_column2->end();
    trimmer_row1->end();
    Fl_Pack* trimmer_row2 = create_row(slicer_width, 25);
    create_process_button(180, "Trim", control_cb, CB_TRIM);
    trimmer_row2->end();

    trimmer_pack->end();
    trimmer_frame->end();

    slicer_panel->end();

    // Rigid shift (xy plane) パネルを作成します。
    const int rigid_shift_width = 600;
    const int rigid_shift_height = 140;
    RolloutPanel* rigid_shift_panel =
        new RolloutPanel(window_width, rigid_shift_width, rigid_shift_height,
            "Rigit shift (xy plane)", false);
    rigid_shift_panel->spacing(10);
    rigid_shift_panel->begin();

    // Rigit shift パネルの 1 列目です。
    const int rigid_shift_column_width = rigid_shift_width / 3 - 20;
    const int rigid_shift_column_height = rigid_shift_height;
    Fl_Pack* rigid_shift_column1 =
        create_column(rigid_shift_column_width - 50, rigid_shift_column_height);

    create_vertical_space(10);

    create_float_input(rigid_shift_column_width - 50, 80, "dx:", &rigid_shift_dx);
    create_float_input(rigid_shift_column_width - 50, 80, "dy:", &rigid_shift_dy);
    create_float_input(rigid_shift_column_width - 50, 80, "dz:", &rigid_shift_dz);

    rigid_shift_column1->end();

    // Rigit shift パネルの 2 列目です。
    Fl_Pack* rigid_shift_column2 =
        create_column(rigid_shift_column_width, rigid_shift_column_height);
    create_vertical_space(10);
    create_float_input(rigid_shift_column_width, 60,
        "Separete at (z):", &rigid_shift_z);
    create_vertical_space(10);
    create_process_button(80, "Shift", control_cb, CB_RIGID_SHIFT);
    create_process_button(80, "Relax", control_cb, CB_RIGID_SHIFT_RELAX);
    rigid_shift_column2->end();

    // Rigit shift パネルの 3 列目です。
    Fl_Pack* rigid_shift_column3 =
        create_column(rigid_shift_column_width, rigid_shift_column_height);
    create_vertical_space(10);
    create_int_input(rigid_shift_column_width, 60,
        "times (x):", &rigid_shift_xtimes);
    create_int_input(rigid_shift_column_width, 60,
        "times (y):", &rigid_shift_ytimes);
    create_vertical_space(10);
    create_process_button(80, "Shift X-Y", control_cb, CB_RIGID_SHIFT_XY);
    rigid_shift_column3->end();

    rigid_shift_panel->end();
    
    // CNT compression (wall mode) パネルを作成します。
    const int cnt_compression_width = 450;
    const int cnt_compression_height = 80;
    RolloutPanel* cnt_compression_panel = new RolloutPanel(
        window_width, cnt_compression_width, cnt_compression_height,
        "CNT compression (wall mode)", true);
    cnt_compression_panel->begin();

    // CNT compression パネルの 1 行目です。
    Fl_Pack* cnt_compression_row1 = create_row(cnt_compression_width, 25, 20);

    create_process_button(120, "Set CNT wall", pointer_cb, CNTWALL_DO_ID);
    create_process_button(120, "Read wall data", pointer_cb, CB_CNTWALL_FB);
    create_process_button(120, "Write wall data", pointer_cb, CNTWALL_WRITE_ID);

    cnt_compression_row1->end();

    cnt_compression_panel->end();

    // Edit パネルを作成します。
    const int edit_width = 120;
    const int edit_height = 80;
    FramedPanel* edit_panel =
        new FramedPanel(window_width, edit_width, edit_height, "Edit", true);
    edit_panel->begin();

    // Edit パネルの 1 行目です。
    Fl_Pack* edit_row1 = create_row(edit_width, 25);
    create_window_open_button(80, "Edit atom", edit_atom_window);
    edit_row1->end();
    edit_panel->end();

    // ボタン行を作成します。
    Fl_Pack* button_row = create_row(window_width, 25, 10);
    create_horizontal_space(100);
    create_process_button(120, "Create", pointer_cb, CREATECONFIG_DO_ID);
    create_process_button(120, "Multiply", pointer_cb, MULTIPLYCELL_DO_ID);
    create_window_close_button(120, "Close");
    button_row->end();

    end_window();
}

int CreateConfigWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 変数の値を初期化します。
            createconfig_mode = 1;
            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            createconfig_mode = 0;
            set_atom_color();
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}
