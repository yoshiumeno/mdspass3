#include "BaseWindow.h"
#include "RolloutPanel.h"
#include "custom_gui_functions.h"

BaseWindow::BaseWindow(const int width, const int height, const char* title)
    : Fl_Window(width + 50, height, title), window_width(width), window_height(height)
{
    //this->resizable(*this); // Dec 2025: commented out
    this->size_range(200, 200, window_width + 50, Fl::h());
    this->begin();

    // スクロール可能にします。
    scroll = new Fl_Scroll(0, 0, window_width + 50, window_height);
    scroll->type(Fl_Scroll::VERTICAL_ALWAYS);
#ifdef __APPLE__
    orig_scroll_cb   = scroll->scrollbar.callback();
    orig_scroll_data = scroll->scrollbar.user_data();
    scroll->scrollbar.callback(scroll_cb, this);
    scroll->scrollbar.when(FL_WHEN_CHANGED);
#endif
    scroll->begin();
    // Dec 2025
    this->resizable(scroll);
    // 水平方向のパックを作成します。
    window_horizontal_pack = create_row(window_width + 50, window_height + 50, 10);
    window_horizontal_pack->begin();
    // 左側に 10 ピクセルのスペースを作成します。
    create_horizontal_space(10);

    // 垂直方向のパックを作成します。
    window_pack = create_column(window_width, window_height);
    // 上部に 10 ピクセルのスペースを作成します。
    create_vertical_space(10);
    window_pack->spacing(10);
    window_pack->begin();
}

void BaseWindow::change_pack_height(int diff)
{
    // Dec 2025: modified for correct resizing（Macでの不具合対応）
    //window_horizontal_pack->size(window_horizontal_pack->w(), window_horizontal_pack->h() + diff);
    int new_h = window_horizontal_pack->h() + diff;
    if (new_h < 0) new_h = 0; // avoiding negative value (just in case)
    window_horizontal_pack->size(window_horizontal_pack->w(), new_h);
    // 中身の高さが変わったことをスクロールに知らせる
    if (scroll) {
        scroll->init_sizes(); // スクロール範囲・内部レイアウトの再計算
        scroll->redraw(); // 描画も更新
    }
    // 既存の redraw も残す（中身の再描画）
    window_pack->redraw();
}

void BaseWindow::begin_window()
{
    this->begin();
    scroll->begin();
    window_horizontal_pack->begin();
    window_pack->begin();
}

void BaseWindow::end_window()
{
    window_pack->end();
    window_horizontal_pack->end();
    scroll->end();
    this->end();
}

#ifdef __APPLE__
void BaseWindow::scroll_cb(Fl_Widget* w, void* data)
{
    auto* self = static_cast<BaseWindow*>(data);
    if (!self) return;
    // まず元々の callback を呼び出して、
    // Fl_Scroll が中身をスクロールする処理をちゃんと動かしてあげる
    if (self->orig_scroll_cb) {
        self->orig_scroll_cb(w, self->orig_scroll_data);
    }
    // その上で、こちらのログを出す (for debug)
    /*
    auto* sb = dynamic_cast<Fl_Scrollbar*>(w);
    int value = sb ? sb->value() : 0;
    fprintf(stderr,
            "[BaseWindow] scroll_cb: scrollbar value=%d (self=%p)\n",
            value, (void*)self);
    fflush(stderr);
    */
    // 3) macOS の描画バグ対策：nudge
    self->nudge_layout();
}
void BaseWindow::nudge_layout()
{
    int X = x(), Y = y(), W = w(), H = h();

    // 1px リサイズ（macOS のレイアウト再評価を強制）
    resize(X, Y, W + 1, H);
    resize(X, Y, W,     H);
}
#endif