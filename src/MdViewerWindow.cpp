#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "MdViewerWindow.h"
#include "GlWindow.h"
#include "custom_gui_functions.h"

MdViewerWindow::MdViewerWindow()
    : Fl_Window(600, 600, "MD viewer")
{
    // ウィンドウのリサイズを可能にします。
    this->resizable(this);

    window_width = 600;
    window_height = 600;

    // OpenGL 描画用のウィンドウを作成します。
    gl_window = new GlWindow(0, 0, window_width, window_height);

    // OpenGL の描画を 30 FPS で更新するためのタイマーを設定します。
    Fl::add_timeout(1.0 / 30.0, gl_window_timer_callback, gl_window);
}

void MdViewerWindow::resize(int X, int Y, int W, int H)
{
    Fl_Window::resize(X, Y, W, H);
    if (gl_window)
    {
        gl_window->resize(0, 0, W, H);
    }
}

MdViewerWindow::~MdViewerWindow()
{
    // タイマーを削除します。
    Fl::remove_timeout(gl_window_timer_callback, gl_window);
}
