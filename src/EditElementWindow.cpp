#include "EditElementWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

EditElementWindow::EditElementWindow() : BaseWindow(180, 300, "Edit element")
{
    // 列の作成を開始します。
    Fl_Pack* edit_element_column1 = new Fl_Pack(0, 0, window_width, 350);
    edit_element_column1->begin();

    // ウィジェットを追加していきます。
    create_vertical_space(10);
    create_int_input(150, 80, "1:", &select_atom[0], 5);
    create_int_input(150, 80, "1-idx:", &select_atom_repidx[0], 5);
    create_int_input(150, 80, "2:", &select_atom[1], 5);
    create_int_input(150, 80, "2-idx:", &select_atom_repidx[1], 5);
    create_int_input(150, 80, "3:", &select_atom[2], 5);
    create_int_input(150, 80, "3-idx:", &select_atom_repidx[2], 5);
    create_int_input(150, 80, "4:", &select_atom[3], 5);
    create_int_input(150, 80, "4-idx:", &select_atom_repidx[3], 5);
    create_vertical_space(10);
    create_process_button(80, "Take this", pointer_cb, EDIT_ELEM_TAKE_ID);
    create_window_close_button(80, "Close");
    edit_element_column1->end();

    end_window();
}

int EditElementWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 値の初期化です。
            for (int i = 0; i < 10; i++) { select_atom[i] = 0; select_atom_repidx[i] = 0; }
            edit_elem_mode = 1;
            draw_replica = 1;

            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            edit_elem_mode = 0;
            draw_replica = 0;
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}

