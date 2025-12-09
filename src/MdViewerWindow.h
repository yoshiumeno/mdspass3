#pragma once

#include <FL/Fl_Window.H>

// GlWindow クラスの前方宣言です。
class GlWindow;

/**
 * MD viewer ウィンドウです。
 **/
class MdViewerWindow : public Fl_Window
{
public:
    /* コンストラクターです。*/
    MdViewerWindow();
    /* デストラクターです。*/
    virtual ~MdViewerWindow();

private:
    /* OpenGL 描画用のウィンドウです。*/
    GlWindow* gl_window;
    /* ウィンドウをリサイズします。*/
    void resize(int X, int Y, int W, int H);
    /* ウィンドウの幅です。*/
    int window_width;
    /* ウィンドウの高さです。*/
    int window_height;
};
