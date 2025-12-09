#pragma once

#include "BaseWindow.h"

/**
 * Edit element ウィンドウです。
 **/
class EditElementWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    EditElementWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
