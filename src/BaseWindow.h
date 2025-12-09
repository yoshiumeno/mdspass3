#pragma once

#include <memory>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Window.H>

/**
 * ウィンドウの基底クラスです。
 * このクラスを継承して、ウィンドウを作成します。
 **/
class BaseWindow : public Fl_Window
{
public:
    /**
     * コンストラクターです。
     *
     * @param width ウィンドウの幅
     * @param height ウィンドウの高さ
     * @param title ウィンドウのタイトル
     **/
    BaseWindow(int width, int height, const char* title);

    /**
     * ウィンドウのパックの高さを変更します。
     *
     * @param diff 高さの変化量
     **/
    void change_pack_height(int diff);

protected:
    /* スクロールです。*/
    Fl_Scroll* scroll;
#ifdef __APPLE__
    Fl_Callback* orig_scroll_cb = nullptr;
    void*        orig_scroll_data = nullptr;
    static void scroll_cb(Fl_Widget* w, void* data);
    void nudge_layout();
#endif
    /* 横方向のパックです。*/
    Fl_Pack* window_horizontal_pack;
    /* 縦方向のパックです。このパックの内部にウィジェットを配置します。*/
    Fl_Pack* window_pack;
    /* ウィンドウの幅です。*/
    int window_width;
    /* ウィンドウの高さです。*/
    int window_height;
    /* ウィンドウの設定を開始します。*/
    void begin_window();
    /* ウィンドウの設定を終了します。*/
    void end_window();
};
