#pragma once

#include <filesystem>
#include "BaseWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

/**
 * Open file ウィンドウです。
 **/
class OpenFileWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    OpenFileWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* ダイアログで表示するディレクトリの文字列を保持する変数です。*/
    std::string _directory_path;
    /* ディレクトリ パスのポインターです。*/
    std::string* _directory_path_ptr = nullptr;
    /**
     * 設定ファイルを読み込みます。
     * @param widget ファイル ブラウザ
     * @param path_ptr ディレクトリ パスのポインター
     **/
    static void read_config_file(Fl_Widget* widget, void* path_ptr);
};
