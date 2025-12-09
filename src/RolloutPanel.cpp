#include "RolloutPanel.h"
#include "BaseWindow.h"

RolloutPanel::RolloutPanel(
    int parent_width,
    int panel_width,
    int panel_height,
    const char* title,
    bool is_vertical,
    bool is_open)
    : FramedPanel(parent_width, panel_width, panel_height, title, is_vertical), _is_open(is_open)
{
    // トグル ボタンにコールバックを追加します。
    _toggle_button->callback(toggle_button_callback, this);
    // トグル ボタンを表示します。
    _toggle_pack->show();

    // タイトル テキストを非表示にします。
    _title_text->hide();

    // パネルの表示、非表示を設定します。
    if (_is_open)
    {
        _horizontal_pack->show();
    }
    else
    {
        _horizontal_pack->hide();
    }
}

void RolloutPanel::toggle()
{
    // パネルの表示、非表示を切り替えます。
    if (_is_open)
    {
        _horizontal_pack->hide();
    }
    else
    {
        _horizontal_pack->show();
    }
    _is_open = !_is_open;

    // 再描画します。
    redraw();
}

void RolloutPanel::redraw()
{
    // ロールアウト パネルの高さの変化量を計算します。
    int diff_height = 0;

    // パネルが開かれた場合です。
    if (_is_open)
    {
        diff_height = _horizontal_pack->h() + _spacing_length;
    }
    // パネルが閉じられた場合です。
    else
    {
        diff_height = -_horizontal_pack->h() - _spacing_length;
    }

    BaseWindow* window = static_cast<BaseWindow*>( _group->window() );
    if (window)
    {
        // ウィンドウのパックの高さを変更します。
        window->change_pack_height(diff_height);
        // ウィンドウを再描画します。
        window->redraw();
    }
}

void toggle_button_callback(Fl_Widget* widget, void* v)
{
    RolloutPanel* panel = static_cast<RolloutPanel*>( v );
    panel->toggle();
}
