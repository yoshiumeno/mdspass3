#include "OpenFileWindow.h"

// 関数のプロトタイプ宣言です。
void readconfig(const char* fname);
void readconfig(const char* fname, int mode);
void readqcelement(const char* fname);
void readqcelement();
void writedata_initialize();

OpenFileWindow::OpenFileWindow() : BaseWindow(600, 500, "Enter filename:")
{
    // 上部に 10 ピクセルのスペースを作成します。
    create_vertical_space(10);

    // 現在のディレクトリを取得します。
    _directory_path = std::filesystem::current_path().string();
    _directory_path_ptr = &_directory_path;
    // ファイル ブラウザを作成します。
    create_file_browser(window_width - 40, window_height - 100, _directory_path, OpenFileWindow::read_config_file, _directory_path_ptr);

    Fl_Pack* row = create_row(600, 25, 20);
    // チェック ボックスを作成します。
    create_check_button(100, "Merge", &imerge, 1);

    // 閉じるボタンを作成します。
    create_window_close_button(100, "Close");

    row->end();
}

int OpenFileWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 何もしません。
            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            // 何もしません。
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}

void OpenFileWindow::read_config_file(Fl_Widget* widget, void* path_ptr)
{
    // 引数をキャストします。
    Fl_File_Browser* browser = static_cast<Fl_File_Browser*>( widget );
    std::string* directory_path_ptr = static_cast<std::string*>( path_ptr );
    std::string directory_path = *directory_path_ptr;

    // ファイル名を保持する配列を作成します。
    char fname[300] = "aaa";
    // 選択されたアイテムのファイル名を取得します。
    std::string fnames;
#if defined(_WIN32) || defined(_WIN64)
    fnames = utf8_to_shift_jis(browser->text(browser->value()));
#else
    fnames = browser->text(browser->value());
#endif
    // 選択されたアイテムのパスを作成します。
    // OpneFileWindow クラスのインスタンスのディレクトリ パスとファイル名を連結します。
    std::filesystem::path selected_path(directory_path + "/" + fnames);
    // 絶対パスに変換します。
    std::filesystem::path absolute_path = std::filesystem::absolute(selected_path);

    strcpy(fname, absolute_path.string().c_str());
    readconfig(fname);
    if (current_qcelement_name != nullptr)
    {
        readqcelement(current_qcelement_name);
    }

    writedata_initialize();
}   
