#pragma once

#include <FL/Fl_Gl_Window.H>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>

/**
 * OpenGL 描画用のウィンドウです。
 **/
class GlWindow : public Fl_Gl_Window
{
public:
    /**
    * コンストラクターです。
    *
    * @param x_ ウィンドウの x 座標
    * @param y_ ウィンドウの y 座標
    * @param w_ ウィンドウの幅
    * @param h_ ウィンドウの高さ
    * @param l ウィンドウのタイトル
    **/
    GlWindow(int x_, int y_, int w_, int h_, const char* l = 0);

    /* ウィンドウをリサイズします。*/
    void resize(int x, int y, int w, int h);

private:
    /* 描画します。*/
    void draw();
    /* ウィンドウ内のイベントを処理します。*/
    int handle(int event);
    /* ウィンドウ内でのクリックを処理します。*/
    void perform_picking(int x, int y, GLuint* selection, GLint hits);
    /* ウィンドウ内でのクリックの終了処理です。*/
    void end_picking(int x, int y, GLuint* selection, GLint hits);
    /* 遠近法あり、なしの描画モードを切り替えます。*/
    void change_ortho(int x, int y, int w, int h);
};
