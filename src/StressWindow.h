#pragma once

#include "BaseWindow.h"
#include "custom_gui_functions.h"

/**
 * Stress ウィンドウです。
 **/
class StressWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    StressWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
