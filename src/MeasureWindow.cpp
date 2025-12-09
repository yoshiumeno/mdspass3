#include "MeasureWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

// 関数のプロトタイプです。
void set_atom_color();

MeasureWindow::MeasureWindow() : BaseWindow(180, 350, "Measure tool")
{
    // 列の作成を開始します。
    Fl_Pack* measure_column1 = new Fl_Pack(0, 0, window_width, 350);
    measure_column1->begin();

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
    create_process_button(80, "Distance 1-2", pointer_cb, MEASURE_CALC2_ID);
    create_process_button(80, "Angle 1-2-3", pointer_cb, MEASURE_CALC3_ID);
    create_process_button(80, "Clear", pointer_cb, MEASURE_CLEAR_ID);
    create_window_close_button(80, "Close");
    measure_column1->end();

    end_window();
}

int MeasureWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 値を初期化します。
            for (int i = 0; i < 10; i++) { select_atom[i] = 0; select_atom_repidx[i] = 0; }
            measure_mode = 1;
            draw_replica = 1;
            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            measure_mode = 0;
            draw_replica = 0;
            set_atom_color();
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}

