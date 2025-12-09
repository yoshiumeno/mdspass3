#pragma once

#include "BaseWindow.h"

/**
 * Edit element ウィンドウです。
 **/
class QuitWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    QuitWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
