#pragma once

#include <string>
#include "BaseWindow.h"

/**
 * Potential file ウィンドウです。
 **/
class PotentialFileWindow :
    public BaseWindow
{
public:
    /* コンストラクターです。*/
    PotentialFileWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* ファイル ブラウザに表示するディレクトリのパスです。*/
    std::string _directory_path;
    /* ディレクトリ パスのポインターです。*/
    std::string* _directory_path_ptr;
    /**
     * pot ファイルを読み込みます。
     * @param widget ファイル ブラウザ
     * @param path_ptr ディレクトリ パスのポインター
     */
    static void read_pot_file(Fl_Widget* widget, void* path_ptr);
};
