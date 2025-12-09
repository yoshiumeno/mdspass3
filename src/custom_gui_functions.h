/* FLTK のウィジェット作成についての関数です。*/
#pragma once

#include <string>
#include <memory>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Spinner.H>
#include <FL/Fl_Radio_Round_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_File_Browser.H>
#include "BaseWindow.h"
#include "CustomChoice.h"
#include "CustomInput.h"

// control_cb のプロトタイプ宣言です。
void control_cb(Fl_Widget* widget, long process_id);

/**
 * パネルを作成します。
 * 親ウィンドウの幅と、パネルの幅を指定し、
 * 親ウィンドウの中央にパネルを配置します。
 *
 * @param parent_width 親ウィンドウの幅
 * @param panel_width パネルの幅
 * @param panel_height パネルの高さ
 * @param is_vertical パネルの方向
 * @return パックのインスタンスへのポインター
 **/
Fl_Pack* create_panel_pack(int parent_width, int panel_width, int panel_height, bool is_vertical);

/**
 * パネルを作成します。
 * x、y 座標を指定し、パネルを配置します。
 *
 * @param x パネルの x 座標
 * @param y パネルの y 座標
 * @param panel_width パネルの幅
 * @param panel_height パネルの高さ
 * @param is_vertical パネルの方向
 * @return パックのインスタンスへのポインター
 **/
Fl_Pack* create_panel_pack_xy(int x, int y, int panel_width, int panel_height, bool is_vertical);

/**
 * 枠付きのグループを作成します。
 * 親ウィンドウの幅と、グループの幅を指定し、
 * 親ウィンドウの中央にグループを配置します。
 *
 * @param parent_width 親ウィンドウの幅
 * @param panel_width パネルの幅
 * @param panel_height パネルの高さ
 * @param title グループのタイトル
 * @return グループのインスタンスへのポインター
 **/
Fl_Group* create_framed_group(int parent_width, int panel_width, int panel_height, const char* title);

/**
 * 枠付きのグループを作成します。
 * x、y 座標を指定し、グループを配置します。
 *
 * @param x グループの x 座標
 * @param y グループの y 座標
 * @param panel_width パネルの幅
 * @param panel_height パネルの高さ
 * @param title グループのタイトル
 * @return グループのインスタンスへのポインター
 **/
Fl_Group* create_framed_group_xy(int x, int y, int panel_width, int panel_height, const char* title);

/**
 * 垂直方向のパックを作成します。
 *
 * @param width パックの幅
 * @param height パックの高さ
 * @param spacing パックの間隔
 * @return パックのインスタンスへのポインター
 **/
Fl_Pack* create_column(int width, int height, int spacing = 0);

/**
 * 水平方向のパックを作成します。
 *
 * @param width パックの幅
 * @param height パックの高さ
 * @param spacing パックの間隔
 * @return パックのインスタンスへのポインター 
 **/
Fl_Pack* create_row(int width, int height, int spacing = 0);

/**
 * 垂直方向のスペースを作成します。
 **/
void create_vertical_space(int height);

/**
 * 水平方向のスペースを作成します。
 **/
void create_horizontal_space(int width);

/**
 * 整数値を入力するスピナーを作成します。
 *
 * @param total_width テキストとスピナーの合計の幅
 * @param spinner_width スピナーの幅
 * @param text テキスト
 * @param int_data int 型の変数へのポインター
 * @param min スピナーの最小値
 * @param max スピナーの最大値
 * @param process_id コールバック関数に渡す処理 ID
 * @return スピナーのインスタンスへのポインター
 **/
Fl_Spinner* create_int_spinner(
    int total_width,
    int spinner_width,
    const char* text,
    int* data,
    int min = 0,
    int max = 10000,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * スピナーのコールバック関数です。
 * 変数へのポインターを受け取り、スピナーの値を変数に反映します。
 *
 * @param widget Fl_Spinner
 * @param data int 型の変数へのポインター
 **/
void int_spinner_callback(Fl_Widget* widget, void* data);

/**
 * 浮動小数点を入力するスピナーを作成します。
 *
 * @param total_width テキストとスピナーの合計の幅
 * @param spinner_width スピナーの幅
 * @param text テキスト
 * @param float_data float 型の変数へのポインター
 * @param min スピナーの最小値
 * @param max スピナーの最大値
 * @param step スピナーのステップ
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return スピナーのインスタンスへのポインター
 **/
Fl_Spinner* create_float_spinner(
    int parent_width,
    int spinner_width,
    const char* text,
    float* data,
    double min = 0,
    double max = 10000,
    double step = 0.5,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * スピナーのコールバック関数です。
 * 変数へのポインターを受け取り、スピナーの値を変数に反映します。
 *
 * @param widget Fl_Spinner
 * @param data float 型の変数へのポインター
 **/
void float_spinner_callback(Fl_Widget* widget, void* data);

/**
 * ドロップダウン リストを作成します。
 *
 * @param total_width ドロップダウン リストの幅
 * @param dropdown_width ドロップダウン リスト内の各項目の幅
 * @param text ラベルテキスト
 * @param selected_index 選択された項目のインデックスを格納するための int 型の変数へのポインター
 * @param item_list ドロップダウン リストに表示する項目のリスト
 * @param item_list_length リストの長さ
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return ドロップダウン リストのインスタンスへのポインター
 **/
CustomChoice* create_dropdown_list(
    int total_width,
    int dropdown_width,
    const char* text,
    int* selected_index,
    const char** item_list = nullptr,
    int item_list_length = -1,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * ドロップダウン リストを作成します。
 * アイテム選択時に potential_set 関数を実行します。
 *
 * @param total_width ドロップダウン リストの幅
 * @param dropdown_width ドロップダウン リスト内の各項目の幅
 * @param text ラベルテキスト
 * @param selected_index 選択された項目のインデックスを格納するための int 型の変数へのポインター
 * @return ドロップダウン リストのインスタンスへのポインター
 **/
CustomChoice* create_dropdown_list_for_potential_set(
    int total_width,
    int dropdown_width,
    const char* text,
    int* selected_index
);

/**
 * ドロップダウン リストのコールバック関数です。
 * 変数へのポインターを受け取り、ドロップダウン リストの選択値を変数に反映します。
 *
 * @param widget Fl_Choice
 * @param data int 型の変数へのポインター
 **/
void dropdown_list_callback(Fl_Widget* widget, void* data);

/**
 * ドロップダウン リストのコールバック関数です。
 * 変数へのポインターを受け取り、ドロップダウン リストの選択値を変数に反映します。
 * potential_set 関数を実行します。
 *
 * @param widget Fl_Choice
 * @param data int 型の変数へのポインター
 **/
void dropdown_list_callback_for_potential_set(Fl_Widget* widget, void* data);

/**
 * チェック ボックスを作成します。
 *
 * @param button_width ボックスとテキストの合計の幅
 * @param text テキスト
 * @param state チェック ボックスと紐づける int 型の変数へのポインター
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return チェック ボックスのインスタンスへのポインター
 **/
Fl_Check_Button* create_check_button(
    int button_width,
    const char* text,
    int* data,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * チェック ボックスのコールバック関数です。
 * 変数へのポインターを受け取り、チェック ボックスの状態を変数に反映します。
 *
 * @param widget Fl_Check_Button
 * @param data int 型の変数へのポインター
 **/
void check_button_callback(Fl_Widget* widget, void* data);

/**
 * テキスト入力フィールドを作成します。
 *
 * @param total_width テキストとテキスト入力フィールドの合計の幅
 * @param input_width テキスト入力フィールドの幅
 * @param label ラベル
 * @param input_text テキスト入力フィールドの値を格納するための char 型の変数へのポインター
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return 入力フィールドのインスタンスへのポインター
 **/
CustomInput* create_text_input(
    int total_width,
    int input_width,
    const char* label,
    char* input_text,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * テキスト入力フィールドのコールバック関数です。
 *
 * @param widget Fl_Input
 * @param data char 型の変数へのポインター
 **/
void text_input_callback(Fl_Widget* widget, void* data);

/**
 * 整数値を入力するテキスト入力フィールドを作成します。
 *
 * @param total_width テキストとテキスト入力フィールドの合計の幅
 * @param input_width テキスト入力フィールドの幅
 * @param label ラベル
 * @param input_value テキスト入力フィールドの値を格納するための int 型の変数へのポインター
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return 入力フィールドのインスタンスへのポインター
 **/
CustomInput* create_int_input(
    int total_width,
    int input_width,
    const char* text,
    int* input_value,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * 浮動小数点を入力するテキスト入力フィールドを作成します。
 *
 * @param total_width テキストとテキスト入力フィールドの合計の幅
 * @param input_width テキスト入力フィールドの幅
 * @param label ラベル
 * @param input_value テキスト入力フィールドの値を格納するための float 型の変数へのポインター
 * @param process_id コールバック関数に渡す処理 ID
 * @return 入力フィールドのインスタンスへのポインター
 **/
CustomInput* create_float_input(
    int total_width,
    int input_width,
    const char* text,
    float* input_value,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * 整数入力フィールドのコールバック関数です。
 *
 * @param widget Fl_Input
 * @param data int 型の変数へのポインター
 **/
void int_input_callback(Fl_Widget* widget, void* data);

/**
 * 浮動小数点入力フィールドのコールバック関数です。
 *
 * @param widget Fl_Input
 * @param data float 型の変数へのポインター
 *
 **/
void float_input_callback(Fl_Widget* widget, void* data);

/**
 * カスタム ラジオ ボタン クラスです。
 * id メンバを持つラジオ ボタンです。
 **/
class CustomRadio : public Fl_Radio_Round_Button
{
public:
    /**
     * コンストラクターです。
     * @param x ボタンの x 座標
     * @param y ボタンの y 座標
     * @param w ボタンの幅
     * @param h ボタンの高さ
     * @param l ボタンのラベル
     * @param id ラジオ ボタンの ID
     **/
    CustomRadio(int x, int y, int w, int h, const char* l, int id)
        : Fl_Radio_Round_Button(x, y, w, h, l), id(id)
    {}

    /* ラジオ ボタンの ID です。*/
    int id;
};

/**
 * ラジオ ボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param text テキスト
 * @param selected_id ラジオ ボタンの状態と紐づける int 型の変数へのポインター
 * @param id ラジオ ボタンの process_id
 * @param process_id コールバック関数に渡す処理 ID
 * @param callback コールバック関数
 * @return ラジオ ボタンのインスタンスへのポインター
 **/
CustomRadio* create_radio_button(
    int button_width,
    const char* text,
    int* data,
    int radio_id,
    int ID = -1,
    Fl_Callback1 callback = control_cb
);

/**
 * ラジオ ボタンのコールバック関数です。
 * 変数へのポインターを受け取り、ラジオ ボタンの状態を変数に反映します。
 *
 * @param widget CustomRadio 型のウィジェット
 * @param data int 型の変数へのポインター
 **/
void radio_button_callback(Fl_Widget* widget, void* data);

/**
 * ボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param label ボタンのラベル
 * @param callback ボタンが押された時に呼び出されるコールバック関数
 * @param process_id コールバック関数に渡す処理 ID
 * @return ボタンのインスタンスへのポインター
 **/
Fl_Button* create_process_button(
    int button_width,
    const char* label,
    Fl_Callback1 callback,
    int process_id = -1
);

/**
 * ウィンドウを開くコールバック関数です。
 *
 * @param widget Fl_Widget 型のインスタンスへのポインター
 * @param new_window void* BaseWindow クラスの派生クラスのインスタンスへのポインター
 */
void open_window_callback(Fl_Widget* widget, void* new_window);

/**
 * ウィンドウを開くボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param label ボタンのラベル
 * @param new_window BaseWindow クラスの派生クラスのインスタンスへのポインター
 **/
template <typename T>
Fl_Button* create_window_open_button(int button_width, const char* label, std::unique_ptr<T>& new_window)
{
    static_assert( std::is_base_of<BaseWindow, T>::value, "T は BaseWindow の派生クラスでなければなりません。" );
    Fl_Button* button = new Fl_Button(0, 0, button_width, 25, label);

    // コールバック関数を設定します。
    button->callback(open_window_callback, new_window.get());

    return button;
}

/**
 * ウィンドウを閉じるボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param label ボタンのラベル
 * @param callback ボタンが押された時に呼び出されるコールバック関数
 * @param process_id コールバック関数に渡す処理 ID
 * @return ボタンのインスタンスへのポインター
 **/
Fl_Button* create_window_close_button(
    int button_width,
    const char* label,
    Fl_Callback1 callback = nullptr,
    int process_id = -1
);

/**
 * ウィンドウを閉じるコールバック関数です。
 *
 * @param widget Fl_Widget 型のインスタンスへのポインター
 * @param data long コールバック関数に渡すデータ
 **/
void close_window_callback(Fl_Widget* widget, long);

/**
 * コールバック関数に渡すデータを格納するための構造体です。
 *
 * @param callback コールバック関数 (引数が Fl_Widget と void*)
 * @param callback コールバック関数 (引数が Fl_Widget と long)
 * @param id コールバック関数に渡す処理 ID
 * @param title タイトル
 * @param int_value_ptr int 型の変数へのポインター
 * @param float_value_ptr float 型の変数へのポインター
 * @param char_value_ptr char 型の変数へのポインター
 * @param path ディレクトリのパス
 * @param path_ptr ディレクトリのパスへのポインター
 **/
struct FlWidgetCallbackData
{
    Fl_Callback* callback0;
    Fl_Callback1* callback;
    int id;
    const char* title;
    int* int_value_ptr;
    float* float_value_ptr;
    char* char_value_ptr;
    std::string path;
    std::string* path_ptr;
};

/**
 * ファイル選択ボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param label ボタンのラベル
 * @param callback ボタンが押された時に呼び出されるコールバック関数
 * @param title ファイル ブラウザのタイトル
 * @param process_id コールバック関数に渡す処理 ID
 * @return ボタンのインスタンスへのポインター
 **/
Fl_Button* create_file_chooser_button(
    int button_width,
    const char* label,
    Fl_Callback1 callback,
    const char* title,
    int process_id = -1
);

/**
 * ファイル選択ボタンのコールバック関数です。
 *
 * @param widget Fl_Widget 型のインスタンスへのポインター
 * @param data void* コールバック関数に渡すデータ
 **/
void file_chooser_button_callback(Fl_Widget* widget, void* data);

/**
 * OpenGL 画面の再描画を行うコールバック関数です。
 * 再描画は 30 FPS で実行します。
 *
 * @param window void* Fl_Gl_Window を継承したクラスのインスタンスへのポインター
 **/
void gl_window_timer_callback(void* window);

/**
 * コールバック関数とその引数を指定したボタンを作成します。
 *
 * @param button_width ボタンの幅
 * @param label ボタンのラベル
 * @param callback ボタンが押された時に呼び出されるコールバック関数
 * @param data コールバック関数に渡す引数
 * @return ボタンのインスタンスへのポインター
 **/
Fl_Button* create_button_with_callback(int button_width, const char* label, Fl_Callback1 callback);

/**
 * キャプチャのコールバック関数です。
 *
 * @param widget Fl_Widget 型のインスタンスへのポインター
 * @param data long コールバック関数に渡すデータ
 **/
void capture_callback(Fl_Widget* widget, long);

/**
 * ファイル ブラウザを作成します。
 *
 * @param width ファイル ブラウザの幅
 * @param height ファイル ブラウザの高さ
 * @param path ファイル ブラウザの初期パス
 * @param callback ファイルが選択された時に呼び出されるコールバック関数
 * @param process_id コールバック関数に渡す処理 ID
 * @return ファイル ブラウザに表示するディレクトリのパス文字列へのポインター
 **/
Fl_File_Browser* create_file_browser(
    int width,
    int height,
    std::string path,
    Fl_Callback callback,
    std::string* path_ptr
);

/**
 * ファイル ブラウザのコールバック関数です。
 *
 * @param widget Fl_Widget 型のインスタンスへのポインター
 * @param data void* コールバック関数に渡すデータ
 **/
void file_browser_callback(Fl_Widget* widget, void* data);

// Windows 環境でのみ使用する関数です。
#if defined _WIN32 || defined _WIN64
/**
 * 文字コードを UTF-8 から Shift_JIS に変換します。
 *
 * @param utf8 UTF-8 文字列
 * @return Shift_JIS 文字列
 **/
std::string utf8_to_shift_jis(const std::string& utf8);

/**
 * 文字コードを Shift_JIS から UTF-8 に変換します。
 *
 * @param Shift_JIS 文字列
 * @return utf8 UTF-8 文字列
 **/
std::string shift_jis_to_utf8(const std::string& shift_jis);
#endif // _WIN32 || _WIN64
