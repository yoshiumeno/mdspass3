#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <filesystem>
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_File_Browser.H>
#include "MdViewerWindow.h"
#include "GlWindow.h"
#include "CustomChoice.h"
#include "CustomInput.h"
#include "custom_gui_functions.h"

// 関数のプロトタイプ宣言です。
void capture();
void potential_set(int control);

Fl_Pack* create_panel_pack(int parent_width, int panel_width, int panel_height, bool is_vertical)
{
    // パネルを作成します。
    Fl_Pack* pack = create_panel_pack_xy(
        ( parent_width - panel_width ) / 2 + 20, 0, panel_width - 40, panel_height, is_vertical);

    return pack;
}

Fl_Pack* create_panel_pack_xy(int x, int y, int panel_width, int panel_height, bool is_vertical)
{
    // パネルを作成します。
    Fl_Pack* pack = new Fl_Pack(x, y, panel_width, panel_height - 30);

    if (is_vertical)
    {
        pack->type(Fl_Pack::VERTICAL);
        pack->begin();
        new Fl_Box(0, 0, 0, 10);
    }
    else
    {
        pack->type(Fl_Pack::HORIZONTAL);
        pack->begin();
        new Fl_Box(0, 0, 10, 0);
    }

    return pack;
}

Fl_Group* create_framed_group(int parent_width, int panel_width, int panel_height, const char* title)
{
    return create_framed_group_xy(
        ( parent_width - panel_width ) / 2, 0, panel_width, panel_height, title);
}

Fl_Group* create_framed_group_xy(int x, int y, int panel_width, int panel_height, const char* title)
{
    // グループを作成します。
    Fl_Group* panel_group = new Fl_Group(x, y, panel_width, panel_height);
    panel_group->begin();

    int frame_y = 30;

    // title が空文字列の場合、frame_y を 0 にします。
    if (std::string(title).empty())
    {
        frame_y = 0;
    }

    // 枠を書きます。
    Fl_Box* frame = new Fl_Box(x, frame_y, panel_width, panel_height - 30, title);
    frame->align(FL_ALIGN_TOP_LEFT);
    frame->box(FL_ENGRAVED_FRAME);
    frame->color(FL_LIGHT2);

    return panel_group;
}

Fl_Pack* create_column(int width, int height, int spacing)
{
    Fl_Pack* pack = new Fl_Pack(0, 0, width, height);
    pack->type(Fl_Pack::VERTICAL);
    pack->spacing(spacing);
    pack->begin();

    return pack;
};

Fl_Pack* create_row(int width, int height, int spacing)
{
    Fl_Pack* pack = new Fl_Pack(0, 0, width, height);
    pack->type(Fl_Pack::HORIZONTAL);
    pack->spacing(spacing);
    pack->begin();

    return pack;
};

void create_vertical_space(int height)
{
    new Fl_Box(0, 0, 0, height);
};

void create_horizontal_space(int width)
{
    new Fl_Box(0, 0, width, 0);
};

void int_spinner_callback(Fl_Widget* widget, void* data)
{

    // 引数をキャストします。
    Fl_Spinner* spinner = static_cast<Fl_Spinner*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // スピナーの値を変数に反映します。
    int* value_ptr = callback_data->int_value_ptr;
    *value_ptr = spinner->value();

    // コールバック関数を実行します。
    callback_data->callback(widget, callback_data->id);
}

void float_spinner_callback(Fl_Widget* widget, void* data)
{

    // 引数をキャストします。
    Fl_Spinner* spinner = static_cast<Fl_Spinner*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // スピナーの値を変数に反映します。
    float* value_ptr = callback_data->float_value_ptr;
    *value_ptr = spinner->value();

    // コールバック関数を実行します。
    callback_data->callback(widget, callback_data->id);
}

Fl_Spinner* create_int_spinner(
    int total_width,
    int spinner_width,
    const char* text,
    int* int_data,
    int min,
    int max,
    int process_id,
    Fl_Callback1 callback)
{
    Fl_Pack* pack = create_row(total_width, 25);

    // テキストを作成します。
    ( new Fl_Box(0, 0, total_width - spinner_width, 25, text) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // スピナーを作成します。
    Fl_Spinner* spinner = new Fl_Spinner(0, 0, spinner_width, 25);

    // スピナーのタイプを整数に設定します。
    spinner->type(FL_INT_INPUT);

    // 変数の値をスピナーの値に反映します。
    spinner->value(*int_data);

    // スピナーの最小値と最大値を設定します。
    spinner->minimum(min);
    spinner->maximum(max);

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->int_value_ptr = int_data;
    data->callback = callback;
    data->id = process_id;

    // コールバック関数を設定します。
    spinner->callback(int_spinner_callback, data);
    pack->end();

    return spinner;
};

Fl_Spinner* create_float_spinner(
    int total_width,
    int spinner_width,
    const char* text,
    float* float_data,
    double min,
    double max,
    double step,
    int process_id,
    Fl_Callback1 callback)
{
    Fl_Pack* pack = create_row(total_width, 25);

    // テキストを作成します。
    ( new Fl_Box(0, 0, total_width - spinner_width, 25, text) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // スピナーを作成します。
    Fl_Spinner* spinner = new Fl_Spinner(0, 0, spinner_width, 25);

    // スピナーのタイプを浮動小数点数に設定します。
    spinner->type(FL_FLOAT_INPUT);

    // 変数の値をスピナーの値に反映します。
    spinner->value(*float_data);

    // スピナーの最小値と最大値を設定します。
    spinner->minimum(min);
    spinner->maximum(max);

    // スピナーのステップを設定します。
    spinner->step(step);

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->float_value_ptr = float_data;
    data->callback = callback;
    data->id = process_id;

    // コールバック関数を設定します。
    spinner->callback(float_spinner_callback, data);
    pack->end();

    return spinner;
};

Fl_Check_Button* create_check_button(
    int button_width,
    const char* text,
    int* state,
    int process_id,
    Fl_Callback1 callback)
{
    Fl_Check_Button* check_button = new Fl_Check_Button(0, 0, button_width, 25, text);

    // 変数の値をチェック ボックスの状態に反映します。
    check_button->value(*state);

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->int_value_ptr = state;
    data->callback = callback;
    data->id = process_id;

    // コールバック関数を設定します。
    check_button->callback(check_button_callback, data);

    return check_button;
}

void check_button_callback(Fl_Widget* widget, void* data)
{

    // 引数をキャストします。
    Fl_Check_Button* check_button = static_cast<Fl_Check_Button*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // チェック ボックスの状態を変数に反映します。
    int* value_ptr = callback_data->int_value_ptr;
    *value_ptr = check_button->value();

    // コールバック関数を実行します。
    if (callback_data->callback)
    {
        callback_data->callback(widget, callback_data->id);
    }
}

CustomRadio* create_radio_button(
    int button_width,
    const char* text,
    int* selected_id,
    int radio_id,
    int process_id,
    Fl_Callback1 callback)
{
    CustomRadio* radio = new CustomRadio(0, 0, button_width, 25, text, radio_id);

    // 変数の値をラジオ ボタンの状態に反映します。
    if (*selected_id == radio_id)
    {
        radio->value(1);
    }

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->int_value_ptr = selected_id;
    data->callback = callback;
    data->id = process_id;

    // コールバック関数を設定します。
    radio->callback(radio_button_callback, data);

    return radio;
};

void radio_button_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    CustomRadio* radio = static_cast<CustomRadio*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // ラジオ ボタンの状態を変数に反映します。
    int* selected_id = callback_data->int_value_ptr;
    if (radio->value())
    {
        *selected_id = radio->id;
    }

    // コールバック関数を実行します。
    if (callback_data->callback)
    {
        callback_data->callback(widget, callback_data->id);
    }
};

CustomChoice* create_dropdown_list(
    int total_width,
    int dropdown_width,
    const char* text,
    int* selected_index,
    const char** item_list,
    int item_list_length,
    int process_id,
    Fl_Callback1 callback)
{
    Fl_Pack* pack = create_row(total_width, 25);

    // テキストを作成します。
    ( new Fl_Box(0, 0, total_width - dropdown_width, 25, text) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // ドロップダウン リストを作成します。
    CustomChoice* choice = new CustomChoice(0, 0, dropdown_width, 25);
    // 変数の値をドロップダウン リストの選択値に反映します。
    choice->attach(selected_index);
    // 選択値の更新を開始します。
    choice->start_polling();

    // アイテム リストが指定されている場合、ドロップダウン リストにアイテムを追加します。
    if (item_list != nullptr)
    {
        for (int i = 0; i < item_list_length; ++i)
        {
            choice->add(item_list[i]);
        }

        // ポインターから初期選択値を設定します。
        if (*selected_index >= 0 && *selected_index < item_list_length)
        {
            choice->value(*selected_index);
        }
    }

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->int_value_ptr = selected_index;
    data->callback = callback;
    data->id = process_id;

    // コールバック関数を設定します。
    choice->callback(dropdown_list_callback, data);
    pack->end();

    return choice;
};

CustomChoice* create_dropdown_list_for_potential_set(
    int total_width,
    int dropdown_width,
    const char* text,
    int* selected_index
)
{
    Fl_Pack* pack = create_row(total_width, 25);

    // テキストを作成します。
    ( new Fl_Box(0, 0, total_width - dropdown_width, 25, text) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // ドロップダウン リストを作成します。
    CustomChoice* choice = new CustomChoice(0, 0, dropdown_width, 25);
    // 変数の値をドロップダウン リストの選択値に反映します。
    choice->attach(selected_index);
    // 選択値の更新を開始します。
    choice->start_polling();

    // コールバック関数を設定します。
    choice->callback(dropdown_list_callback_for_potential_set, selected_index);

    pack->end();

    return choice;
};

void dropdown_list_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    Fl_Choice* choice = static_cast<Fl_Choice*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // ドロップダウン リストの選択値を変数に反映します。
    int* selected_index = callback_data->int_value_ptr;
    *selected_index = choice->value();

    // コールバック関数を実行します。
    if (callback_data->callback)
    {
        callback_data->callback(widget, callback_data->id);
    }
};

void dropdown_list_callback_for_potential_set(Fl_Widget* widget, void* data)
{
    Fl_Choice* choice = static_cast<Fl_Choice*>( widget );
    int* selected_index = static_cast<int*>( data );
    *selected_index = choice->value();
    // potential_set 関数を実行します。
    potential_set(0);
}

CustomInput* create_text_input(
    int total_width,
    int input_width,
    const char* label,
    char* input_text,
    int process_id,
    Fl_Callback1 callback)
{
    Fl_Pack* pack = create_row(total_width, 25);

    // テキストを作成します。
    ( new Fl_Box(0, 0, total_width - input_width, 25, label) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // テキスト入力フィールドを作成します。
    CustomInput* input = new CustomInput(0, 0, input_width, 25);
    // input_text のサイズを取得します。
    size_t input_text_size = sizeof(input_text);
    // 変数の値をテキスト入力フィールドの値に反映します。
    input->attach(input_text, input_text_size);
    // 入力フィールドの値を更新するためのポーリングを開始します。
    input->start_polling();

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->char_value_ptr = input_text;
    data->callback = callback;
    data->id = process_id;
    void* callback_arg = static_cast<void*>( data );

    // コールバック関数を設定します。
    input->callback(text_input_callback, data);
    pack->end();

    return input;
}

void text_input_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // コールバック関数を実行します。
    if (callback_data->callback)
    {
        callback_data->callback(widget, callback_data->id);
    }
};

CustomInput* create_int_input(
    int total_width,
    int input_width,
    const char* label,
    int* input_value,
    int process_id,
    Fl_Callback1 callback)
{
    // ラベルとテキスト入力フィールドを含むパックを作成します。
    Fl_Pack* pack = create_row(total_width, 25);

    // ラベルを作成します。
    ( new Fl_Box(0, 0, total_width - input_width, 25, label) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // テキスト入力フィールドを作成します。
    CustomInput* input = new CustomInput(0, 0, input_width, 25);
    input->type(FL_INT_INPUT);
    // 変数の値をテキスト入力フィールドの値に反映します。
    input->attach(input_value);
    // 入力フィールドの値を更新するためのポーリングを開始します。
    input->start_polling();

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->int_value_ptr = input_value;
    data->callback = callback;
    data->id = process_id;
    void* callback_arg = static_cast<void*>( data );

    input->callback(int_input_callback, callback_arg);

    pack->end();
    return input;
}

CustomInput* create_float_input(
    int total_width,
    int input_width,
    const char* label,
    float* input_value,
    int process_id,
    Fl_Callback1 callback)
{
    // ラベルとテキスト入力フィールドを含むパックを作成します。
    Fl_Pack* pack = create_row(total_width, 25);

    // ラベルを作成します。
    ( new Fl_Box(0, 0, total_width - input_width, 25, label) )->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    // テキスト入力フィールドを作成します。
    CustomInput* input = new CustomInput(0, 0, input_width, 25);
    // テキスト入力フィールドのタイプを浮動小数点数に設定します。
    input->type(FL_FLOAT_INPUT);
    // 変数の値をテキスト入力フィールドの値に反映します。
    input->attach(input_value);
    // 入力フィールドの値を更新するためのポーリングを開始します。
    input->start_polling();

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->float_value_ptr = input_value;
    data->callback = callback;
    data->id = process_id;
    void* callback_arg = static_cast<void*>( data );

    // コールバック関数を設定します。
    input->callback(float_input_callback, callback_arg);

    pack->end();
    return input;
}

void int_input_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // コールバック関数を実行します。
    if (callback_data->callback)
    {
        callback_data->callback(widget, callback_data->id);
    }
}

void float_input_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );

    // コールバック関数を実行します。
    callback_data->callback(widget, callback_data->id);
}

Fl_Button* create_process_button(int button_width, const char* label, Fl_Callback1 callback, int process_id)
{
    Fl_Button* button = new Fl_Button(0, 0, button_width, 25, label);

    button->argument(process_id);
    long temp_long = button->argument();

    // コールバック関数と処理 ID を設定します。
    button->callback(callback, process_id);

    return button;
}

Fl_Button* create_window_close_button(int button_width, const char* label, Fl_Callback1 callback, int process_id)
{
    Fl_Button* button = new Fl_Button(0, 0, button_width, 25, label);

    // コールバック関数が引数で渡されている場合、それを設定します。
    if (callback)
    {
        button->callback(callback, process_id);
    }
    // コールバック関数が引数で渡されていない場合、デフォルトのコールバック関数を設定します。
    else
    {
        button->callback(close_window_callback, -1);
    }
    return button;
}

void close_window_callback(Fl_Widget* widget, long)
{
    Fl_Window* window = widget->window();
    if (window) window->hide();
}

void open_window_callback(Fl_Widget* widget, void* new_window)
{
    BaseWindow* window = static_cast<BaseWindow*>( new_window );
    if (window) window->show();
}

Fl_Button* create_button_with_callback(int button_width, const char* label, Fl_Callback1 callback)
{
    Fl_Button* button = new Fl_Button(0, 0, button_width, 25, label);

    // コールバック関数と処理 ID を設定します。
    button->callback(callback);

    return button;
}

void capture_callback(Fl_Widget* widget, long)
{
#ifndef NOPNG
    // キャプチャを実行します。
    capture();
#endif
}

Fl_Button* create_file_chooser_button(
    int button_width,
    const char* label,
    Fl_Callback1 callback,
    const char* title,
    int process_id)
{
    Fl_Button* file_chooser_button = new Fl_Button(0, 0, button_width, 25, label);

    // callback と process_id を data にまとめて、ボタンのコールバックの引数にします。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->callback = callback;
    data->id = process_id;
    data->title = title;

    file_chooser_button->callback(file_chooser_button_callback, (void*)data);
    return file_chooser_button;
}

void file_chooser_button_callback(Fl_Widget* widget, void* data)
{
    // ファイル ブラウザに関するデータを取得します。
    FlWidgetCallbackData* file_chooser_data = static_cast<FlWidgetCallbackData*>( data );

    // ファイル選択ダイアログを表示します。
    Fl_File_Chooser* file_chooser = new Fl_File_Chooser(
        ".", "*.*", Fl_File_Chooser::SINGLE, file_chooser_data->title);
    file_chooser->show();
    while (file_chooser->shown())
    {
        Fl::wait();
    }

    // ファイルが選択された場合、コールバック関数を実行します。
    if (file_chooser->value() != NULL)
    {
        if (file_chooser_data->callback)
        {
            file_chooser_data->callback((Fl_Widget*)file_chooser, file_chooser_data->id);
        }
    }

    delete file_chooser;
}

Fl_File_Browser* create_file_browser(int width, int height, std::string path, Fl_Callback callback, std::string* path_ptr)
{
    // ファイル ブラウザを作成します。
    Fl_File_Browser* file_browser = new Fl_File_Browser(0, 0, width, height);
    // ブラウザのタイプを設定します。
    file_browser->type(FL_HOLD_BROWSER);
    // ディレクトリを指定します。
    file_browser->load(".");

    // コールバックの引数となるデータを作成します。
    FlWidgetCallbackData* data = new FlWidgetCallbackData();
    data->callback0 = callback;
    data->path = path;
    data->path_ptr = path_ptr;

    // コールバック関数を設定します。
    file_browser->callback(file_browser_callback, (void*)data);

    return file_browser;
}

void file_browser_callback(Fl_Widget* widget, void* data)
{
    // 引数をキャストします。
    Fl_File_Browser* file_browser = static_cast<Fl_File_Browser*>( widget );
    FlWidgetCallbackData* callback_data = static_cast<FlWidgetCallbackData*>( data );
    std::string* directory_path_ptr = static_cast<std::string*>(callback_data->path_ptr);

    // 選択されたアイテムのパスを取得します。
    std::string selected;
#if defined(_WIN32) || defined(_WIN64)
    const char* selected_windows = file_browser->text(file_browser->value());
    // 未選択の場合、何もしません。
    if (!selected_windows) return;
    selected = utf8_to_shift_jis(selected_windows);
#else
    const char* selected_linux = file_browser->text(file_browser->value());
    // 未選択の場合、何もしません。
    if (!selected_linux) return;
    selected = std::string(selected_linux);
#endif

    // ブラウザが表示しているパスを取得します。
    std::string current_path = callback_data->path;
    // 選択されたアイテムのパスを取得します。
    std::filesystem::path selected_path(current_path + "/" + selected);
    // 絶対パスを取得します。
    std::filesystem::path absolute_path = std::filesystem::absolute(selected_path);

    // ディレクトリを選択した場合、そのディレクトリに移動して、ファイル ブラウザを更新します。
    if (std::filesystem::is_directory(absolute_path))
    {
#if defined(_WIN32) || defined(_WIN64)
        // Windows の場合、文字コードを変換します。
        std::string new_path = shift_jis_to_utf8(absolute_path.string().c_str());
        // ファイル ブラウザに表示するディレクトリを変更します。
        file_browser->load(new_path.c_str());
        // パスを更新します。
        callback_data->path = utf8_to_shift_jis(new_path);
#else
        std::string new_path = absolute_path.string();
        file_browser->load(new_path.c_str());
        callback_data->path = new_path;
#endif
        return;
    }

    // コールバック関数を実行します。
    if (callback_data->callback0)
    {
        // ファイルが選択されている場合のみ、コールバック関数を実行します。
        if (file_browser->value() > 0)
        {
            directory_path_ptr = &callback_data->path;
            callback_data->callback0(widget, directory_path_ptr);
        }
    }
}

void gl_window_timer_callback(void* window)
{
    GlWindow* gl_window = static_cast<GlWindow*>( window );
    // OpenGL ウィンドウを再描画します。
    gl_window->redraw();
    gl_window->resize(gl_window->x(), gl_window->y(), gl_window->w(), gl_window->h());
    // 30 FPS で画面を更新します。
    Fl::repeat_timeout(1.0 / 30.0, gl_window_timer_callback, window);
}

// Windows 環境でのみ以下の関数を定義します。
#if defined(_WIN32) || defined(_WIN64)
std::string utf8_to_shift_jis(const std::string& utf8)
{
    // バッファー サイズを取得します。
    int length_wide = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, nullptr, 0);
    // 取得に失敗した場合、空文字列を返します。
    if (length_wide == 0)
    {
        return "";
    }

    // UTF-16 の文字列を格納するためのバッファを作成します。
    std::wstring utf16(length_wide - 1, L'\0');
    // UTF-8 から UTF-16 へ変換します。
    MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, &utf16[0], length_wide);

    // バッファー サイズを取得します。
    int length_shift_jis = WideCharToMultiByte(CP_ACP, 0, utf16.c_str(), -1, nullptr, 0, nullptr, nullptr);
    // 取得に失敗した場合、空文字列を返します。
    if (length_shift_jis == 0)
    {
        return "";
    }

    // Shift-JIS の文字列を格納するためのバッファを作成します。
    std::string shift_jis(length_shift_jis - 1, '\0');
    // UTF-16 から Shift-JIS へ変換します。
    WideCharToMultiByte(CP_ACP, 0, utf16.c_str(), -1, &shift_jis[0], length_shift_jis, nullptr, nullptr);

    return shift_jis;
}

std::string shift_jis_to_utf8(const std::string& shift_jis)
{
    // バッファー サイズを取得します。
    int length_wide = MultiByteToWideChar(CP_ACP, 0, shift_jis.c_str(), -1, nullptr, 0);
    // 取得に失敗した場合、空文字列を返します。
    if (length_wide == 0)
    {
        return "";
    }

    // UTF-16 の文字列を格納するためのバッファを作成します。
    std::wstring utf16(length_wide - 1, L'\0');
    // Shift-JIS から UTF-16 へ変換を実行します。
    MultiByteToWideChar(CP_ACP, 0, shift_jis.c_str(), -1, &utf16[0], length_wide);

    // バッファー サイズを取得します。
    int length_utf8 = WideCharToMultiByte(CP_UTF8, 0, utf16.c_str(), -1, nullptr, 0, nullptr, nullptr);
    // 取得に失敗した場合、空文字列を返します。
    if (length_utf8 == 0)
    {
        return "";
    }

    // UTF-8 の文字列を格納するためのバッファを作成します。
    std::string utf8(length_utf8 - 1, '\0');
    // UTF-16 から UTF-8 へ変換します。
    WideCharToMultiByte(CP_UTF8, 0, utf16.c_str(), -1, &utf8[0], length_utf8, nullptr, nullptr);

    return utf8;
}
#endif
