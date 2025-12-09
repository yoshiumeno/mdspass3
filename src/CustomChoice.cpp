#include "CustomChoice.h"
#include "myheader.h"

CustomChoice::CustomChoice(int X, int Y, int W, int H, const char* L)
    : Fl_Choice(X, Y, W, H, L), _polling(false), _external_selected_index(nullptr)
{}

CustomChoice::~CustomChoice()
{
    // 値の更新を終了します。
    stop_polling();
}

void CustomChoice::attach(int* selected_index)
{
    // 外部変数を設定します。
    //_external_choices = new_choices;
    _external_selected_index = selected_index;
    // 初期同期を行います。
    sync_with_external_variables();
}

void CustomChoice::start_polling()
{
    if (!_polling && _external_selected_index)
    {
        // 値の更新を開始します。
        Fl::add_timeout(POLLING_INTERVAL, polling_callback, this);
        _polling = true;
    }
}

void CustomChoice::stop_polling()
{
    if (_polling)
    {
        // 値の更新を終了します。
        Fl::remove_timeout(polling_callback, this);
        _polling = false;
    }
}

void CustomChoice::polling_callback(void* v)
{
    // 引数を元の型にキャストします。
    CustomChoice* widget = static_cast<CustomChoice*>( v );
    // 外部変数の値を同期します。
    widget->sync_with_external_variables();
    if (widget->_polling)
    {
        // 一定時間経過後に再度更新します。
        Fl::repeat_timeout(POLLING_INTERVAL, polling_callback, widget);
    }
}

void CustomChoice::sync_with_external_variables()
{
    if (_external_selected_index != nullptr)
    {
        // 選択されたインデックスが選択肢リストの範囲内かチェックします。
        if (*_external_selected_index < static_cast<int>( this->size() ))
        {
            // 新しいインデックスを設定します。
            this->value(*_external_selected_index);
        }
        // ウィジェットを再描画します。
        this->redraw();
    }
}

int CustomChoice::handle(int event)
{
    switch (event)
    {
        case FL_FOCUS:
            // 入力開始時 (フォーカス取得時) に自動更新を停止します。
            stop_polling();
            break;
        case FL_UNFOCUS:
            // 入力終了時 (フォーカス失った時) に変数を更新し、自動更新を再開します。
            sync_with_external_variables();
            start_polling();
            break;
        default:
            break;
    }
    // FLTK のデフォルトのイベント ハンドリングを継続します。
    return Fl_Choice::handle(event);
}
