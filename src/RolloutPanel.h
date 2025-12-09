#pragma once

#include <FL/Fl_Group.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Button.H>
#include "FramedPanel.h"
#include "custom_gui_functions.h"

/**
 * ロールアウト パネルを表すクラスです。
 **/
class RolloutPanel : public FramedPanel
{
public:
    /**
    * コンストラクターです。
    *
    * @param parent_width 親ウィンドウの幅
    * @param width パネルの幅
    * @param height パネルの高さ
    * @param title タイトル
    * @param is_vertical パネル内に配置するウィジェットの配置方向
    * @param is_open パネルを開いた状態にするかどうか
    **/
    RolloutPanel(
        int parent_width,
        int width,
        int height,
        const char* title,
        bool is_vertical,
        bool is_open = false);

    /* トグルします。*/
    void toggle();
    /* 再描画します。*/
    void redraw();

private:
    /* ロールアウトが開かれているかどうかのフラグです。*/
    bool _is_open;
};

/* トグル ボタンのコールバック関数です。*/
void toggle_button_callback(Fl_Widget* widget, void* data);
