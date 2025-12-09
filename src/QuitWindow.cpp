#include "QuitWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

QuitWindow::QuitWindow() : BaseWindow(300, 100, "Quit confirmation")
{
    // 列の作成を開始します。
    Fl_Pack* quit_column1 = new Fl_Pack(0, 0, window_width, 50);
    quit_column1->begin();

    // ウィジェットを追加していきます。
    create_vertical_space(20);
    new Fl_Box(0, 0, 10, 30, "Do you really want to quit?");
    create_vertical_space(20);
    Fl_Pack* quit_row1 = create_row(window_width, 50, 30);
    create_horizontal_space(1);
    create_process_button(110, "Yes", control_cb, CB_QUIT);
    create_window_close_button(110, "No");
    quit_row1->end();
    quit_column1->end();

    end_window();
}

int QuitWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 値の初期化です。
            //for (int i = 0; i < 10; i++) { select_atom[i] = 0; select_atom_repidx[i] = 0; }
            //edit_elem_mode = 1;
            //draw_replica = 1;

            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            //edit_elem_mode = 0;
            //draw_replica = 0;
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}

