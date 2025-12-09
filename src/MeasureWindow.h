#pragma once

#include "BaseWindow.h"

/**
 * Measure ウィンドウです。
 **/
class MeasureWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    MeasureWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
