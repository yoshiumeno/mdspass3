#pragma once

#include <string>
#include <FL/Fl_Group.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Button.H>
#include <Fl/Fl_Box.H>

/**
 * フレーム付きのパネルを表すクラスです。
 **/
class FramedPanel
{
public:
    /**
    * コンストラクターです。
    *
    * @param parent_width 親の幅
    * @param panel_width パネルの幅
    * @param panel_height パネルの高さ
    * @param title パネルのタイトル
    * @param is_vertical 縦方向のパネルかどうか
    **/
    FramedPanel(int parent_width, int panel_width, int panel_height, const char* title, bool is_vertical);
    /* デストラクターです。*/
    virtual ~FramedPanel();
    /* パネルを開始します。*/
    void begin();
    /* パネルを終了します。*/
    void end();
    /* パネル内のオブジェクト間のスペースを指定します。*/
    void spacing(int spacing);

protected:
    /* FramedPanel のグループです。*/
    Fl_Group* _group;
    /* 水平方向のパックです。*/
    Fl_Pack* _horizontal_pack;
    /* 垂直方向のパックです。*/
    Fl_Pack* _pack;
    /* トグル ボタンのパックです。*/
    Fl_Pack* _toggle_pack;
    /* トグル ボタンです。*/
    Fl_Button* _toggle_button;
    /* タイトルです。*/
    Fl_Box* _title_text;
    /* 親の幅です。*/
    int _parent_width;
    /* パネルの幅です。*/
    int _panel_width;
    /* パネルの高さです。*/
    int _panel_height;
    /* パネルのタイトルです。*/
    const char* _title;
    /* パックが縦方向かどうかです。*/
    bool _is_vertical;
    /* スペーシングの値です。*/
    int _spacing_length;
};
