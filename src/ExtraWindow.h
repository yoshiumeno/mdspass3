#pragma once

#include <memory>
#include "BaseWindow.h"
#include "MeasureWindow.h"

/**
 * Extra ウィンドウです。
 **/
class ExtraWindow : public BaseWindow
{
public:
    /*コンストラクターです。*/
    ExtraWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* Measure ウィンドウです。*/
    std::unique_ptr<MeasureWindow> measure_window;
};
