#pragma once

#include <memory>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Pack.H>
#include "BaseWindow.h"
#include "PotentialFileWindow.h"

/**
 * Potential arg ウィンドウです。
 **/
class PotentialArgWindow :
    public BaseWindow
{
public:
    /* コンストラクターです。*/
    PotentialArgWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* Potential file ウィンドウです。*/
    std::unique_ptr<PotentialFileWindow> potential_file_window;
    /* Potential 選択が可能な場合に表示するパックです。*/
    Fl_Pack* _potential_option_pack_select;
    /* Potential 選択ができない場合に表示するパックです。*/
    Fl_Pack* _potential_option_pack_none;
    /* 選択可能な Potential の数です。*/
    int _narg = 0;
    // Potential 選択のドロップダウン リストです。
    Fl_Choice* _arg_choice = nullptr;
    // Potential ファイルが読み込み可能かどうかです。
    bool _potarg_readable = false;
    // Potential ファイル選択が可能な場合に表示するパックです。
    Fl_Pack* _potential_file_pack_select;
    // Potential ファイル選択ができない場合に表示するパックです。
    Fl_Pack* _potential_file_pack_none;
};
