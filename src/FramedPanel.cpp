#include "FramedPanel.h"
#include "custom_gui_functions.h"

FramedPanel::FramedPanel(
    int parent_width,
    int panel_width,
    int panel_height,
    const char* title,
    bool is_vertical)
    :
    _parent_width(parent_width),
    _panel_width(panel_width),
    _panel_height(panel_height),
    _title(title),
    _is_vertical(is_vertical),
    _spacing_length(0)
{
    // トグル ボタンを作成します。
    _toggle_pack = create_row(_parent_width, 25);
    create_horizontal_space(( _parent_width - 250 ) / 2);
    _toggle_button = new Fl_Button(0, 0, 250, 25, _title);
    // トグル ボタンの外観を設定します。
    _toggle_button->box(FL_ROUND_UP_BOX);
    _toggle_pack->end();
    _toggle_pack->hide();

    // タイトルを表示します。
    _title_text = new Fl_Box(0, 0, 280, 25, _title);
    _title_text->align(FL_ALIGN_INSIDE | FL_ALIGN_TOP);

    // パックを作成します。
    _horizontal_pack = create_row(_parent_width, _panel_height);
    create_horizontal_space(( _parent_width - panel_width ) / 2);
    // 枠付きのグループを作成します。
    _group = create_framed_group(_parent_width, _panel_width, _panel_height, "");
    _pack = create_panel_pack(_parent_width, _panel_width, _panel_height, _is_vertical);
}

FramedPanel::~FramedPanel()
{
    delete _toggle_pack;
    delete _horizontal_pack;
}

void FramedPanel::spacing(int spacing)
{
    _pack->spacing(spacing);
    _spacing_length = spacing;
}

void FramedPanel::end()
{
    _pack->end();
    _group->end();
    _horizontal_pack->end();
}

void FramedPanel::begin()
{
    _horizontal_pack->begin();
    _group->begin();
    _pack->begin();
}
