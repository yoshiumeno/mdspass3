#include "PotentialArgWindow.h"
#include "FramedPanel.h"
#include "custom_gui_functions.h"
#include "myheader.h"

PotentialArgWindow::PotentialArgWindow()
    : BaseWindow(350, 400, "Potential arg")
{
    end_window();
    // 別のウィンドウのインスタンスを作成します。
    potential_file_window = std::make_unique<PotentialFileWindow>();
    begin_window();

    // Potential option パネルを作成します。
    const int potential_option_width = window_width - 50;
    const int potential_option_height = 120;

    ::Fl_Group* potential_option_frame = create_framed_group_xy(
        0, 0, potential_option_width, potential_option_height, "Potential argument");

    // Potential を引数から選択する場合に表示するパックです。
    _potential_option_pack_select = create_panel_pack_xy(
        10, 30, potential_option_width, potential_option_height, true);
    create_vertical_space(10);
    new Fl_Box(0, 0, 300, 0, "Select potential argument from list.");
    create_vertical_space(20);
    Fl_Pack* option_select_row1 = create_row(potential_option_width, 25);
    create_horizontal_space((potential_option_width - 200) / 2);
    // Potential argument を選択するドロップダウン リストを作成します。
    _arg_choice = create_dropdown_list_for_potential_set(200, 150, "arg:", &ipotarg);
    option_select_row1->end();
    _potential_option_pack_select->end();
    // パックを非表示にします。
    _potential_option_pack_select->hide();

    // Potential を引数から選択できない場合に表示するパックです。
    _potential_option_pack_none = create_panel_pack_xy(
        10, 30, potential_option_width, potential_option_height, true);
    create_vertical_space(10);
    new Fl_Box(0, 0, 300, 0, "No argument available for this potential!");
    _potential_option_pack_none->end();

    potential_option_frame->end();

    // Potential file パネルを作成します。
    const int potential_file_width = window_width - 50;
    const int potential_file_height = 100;

    ::Fl_Group* potential_file_frame = create_framed_group_xy(
        0, 0, potential_file_width, potential_file_height, "Potential file");

    // Pot ファイル選択ボタンを含むパネルを作成します。
    _potential_file_pack_select = create_panel_pack_xy(
        10, 30, potential_file_width, potential_file_height, true);
    create_vertical_space(10);
    Fl_Pack* row1 = create_row(potential_file_width, 25);
    create_horizontal_space(( potential_file_width - 120 ) / 2);
    create_window_open_button(120, "Pot file", potential_file_window);
    row1->end();
    _potential_file_pack_select->end();
    // パックを非表示にします。
    _potential_file_pack_select->hide();

    // Pot ファイルを読み込めないパラメーターに対して表示するパネルを作成します。
    _potential_file_pack_none = create_panel_pack_xy(
        10, 30, potential_file_width, potential_file_height, true);
    create_vertical_space(10);
    new Fl_Box(0, 0, 300, 0, "No readable file of parameters!");
    _potential_file_pack_none->end();

    potential_file_frame->end();

    // ウィンドウを閉じるボタンを作成します。
    create_window_close_button(120, "Close");

    end_window();
}

int PotentialArgWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // Potential option パネルの表示を制御します。
            _narg = atom.potarg_number[ipottype];
            // Potential argumennt が選択可能な場合です。
            if (_narg > 0 && _arg_choice != nullptr)
            {
                // ドロップダウン リストをクリアします。
                _arg_choice->clear();
                // ドロップダウン リストにアイテムを追加します。
                for (int i = 0; i < _narg; ++i)
                {
                    _arg_choice->add(atom.potargstring_list[ipottype][i]);
                }
                // ドロップダウン リストを含むパックを表示します。
                _potential_option_pack_select->show();
                _potential_option_pack_none->hide();
            }
            // Potential argumennt が選択できない場合です。
            else
            {
                // ドロップダウン リストを含むパックを非表示にします。
                _potential_option_pack_select->hide();
                _potential_option_pack_none->show();
            }

            // Potential file パネルの表示を制御します。
            _potarg_readable = atom.potarg_readable[ipottype];
            // Potential file が読み込める場合です。
            if (_potarg_readable)
            {
                _potential_file_pack_select->show();
                _potential_file_pack_none->hide();
            }
            // Potential file が読み込めない場合です。
            else
            {
                _potential_file_pack_select->hide();
                _potential_file_pack_none->show();
            }
            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            // 何もしません。
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}
