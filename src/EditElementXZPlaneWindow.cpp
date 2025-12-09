#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <string>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Window.H>
#include "EditElementXZPlaneWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

extern GLfloat** color;
extern GLfloat yellow[];

EditElementXZPlaneWindow::EditElementXZPlaneWindow()
    : BaseWindow(180, 300, "Edit element (XZ-Plane Mode)")
{
    // 列の作成を開始します。
    Fl_Pack* edit_element_xz_column1 = create_column(window_width, 350);

    create_vertical_space(10);
    // ウィジェットを追加していきます。
    create_int_input(150, 80, "1:", &select_atom[0], 5);
    create_int_input(150, 80, "1-idx:", &select_atom_repidx[0], 5);
    create_int_input(150, 80, "2:", &select_atom[1], 5);
    create_int_input(150, 80, "2-idx:", &select_atom_repidx[1], 5);
    create_int_input(150, 80, "3:", &select_atom[2], 5);
    create_int_input(150, 80, "3-idx:", &select_atom_repidx[2], 5);
    create_int_input(150, 80, "4:", &select_atom[3], 5);
    create_int_input(150, 80, "4-idx:", &select_atom_repidx[3], 5);
    create_vertical_space(10);
    create_process_button(80, "Take this", pointer_cb, EDIT_ELEM_XZ_TAKE_ID);
    create_window_close_button(80, "Close");
    edit_element_xz_column1->end();

    end_window();
}

int EditElementXZPlaneWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 値の初期化です。
            for (int i = 0; i < 10; i++)
            {
                select_atom[i] = 0;
                select_atom_repidx[i] = 0;
            }
            edit_elem_mode = 1;
            draw_replica = 2;

            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            edit_elem_mode = 0;
            draw_replica = 0;
            for (int i = 1; i <= atom.natom + icnt; i++) { memcpy(::color[i], yellow, sizeof(GLfloat) * 4); }
            return 1;
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}
