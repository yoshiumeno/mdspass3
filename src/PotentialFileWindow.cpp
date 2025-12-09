#include <filesystem>
#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif
#include "PotentialFileWindow.h"
#include "custom_gui_functions.h"
#include "myheader.h"

extern std::string potfile;
// 関数のプロトタイプ宣言です。
void set_potfile();

PotentialFileWindow::PotentialFileWindow() : BaseWindow(600, 500, "Potential parameter filename")
{
    // スペースを追加します。
    create_vertical_space(20);

    // ドロップダウン リストを作成します。
    Fl_Choice* potential_choice = create_dropdown_list_for_potential_set(250, 100, "Potencial:", &ipottype);
    for (int i = 0; i < MAXPOTTYPE; ++i)
    {
        potential_choice->add(atom.potstring_list[i]);
    }
    potential_choice->value(ipottype);

    // スペースを追加します。
    create_vertical_space(20);

    // 現在のディレクトリを取得します。
    _directory_path = std::filesystem::current_path().string();
    _directory_path_ptr = &_directory_path;
    // ファイル ブラウザを作成します。
    create_file_browser(window_width - 40, window_height - 200, _directory_path, PotentialFileWindow::read_pot_file, _directory_path_ptr);

    // 閉じるボタンを作成します。
    create_window_close_button(100, "Close");
}

int PotentialFileWindow::handle(int event)
{
    switch (event)
    {
        // ウィンドウが表示された時の処理です。
        case FL_SHOW:
            // 2024/02 CJS main.cpp からそのまま移動した処理です。
#ifdef __GNUC__
            chdir("pot");
#endif
#ifdef _WIN32
            //SetCurrentDirectory("pot");
            chdir("pot");
#endif
            // 2024/02 CJS ここまで

            return 1;

        // ウィンドウが閉じられた時の処理です。
        case FL_HIDE:
            // 2024/02 CJS main.cpp からそのまま移動した処理です。
#ifdef __GNUC__
            chdir(cwdname);
#endif
#ifdef _WIN32
            //SetCurrentDirectory(cwdname);
            chdir(cwdname);
#endif
            // 2024/02 CJS ここまで
            return 1;

        default:
            return Fl_Window::handle(event);
    }
}

void PotentialFileWindow::read_pot_file(Fl_Widget* widget, void* path_ptr)
{
    // ウィジェットをファイル ブラウザにキャストします。
    Fl_File_Browser* file_browser = (Fl_File_Browser*)( widget );
    std::string* directory_path_ptr = static_cast<std::string*>(path_ptr);
    std::string directory_path = *directory_path_ptr;

    // ファイル ブラウザの選択されたファイル名を取得します。
    potfile += file_browser->text(file_browser->value());
    std::string argtext = file_browser->text(file_browser->value());
    // 選択されたアイテムのパスを作成します。
    std::filesystem::path selected_path(directory_path + "/" + argtext);
    // 絶対パスに変換します。
    std::filesystem::path absolute_path = std::filesystem::absolute(selected_path);
    // ファイル名を Atom クラスの potential_arg メンバにコピーします。
    strcpy(atom.potential_arg, argtext.c_str());
    // pot ファイルを更新します。
    potfile = absolute_path.string();

    set_potfile();
}
