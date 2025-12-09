#pragma once

#include "BaseWindow.h"

/**
* Edit atom ウィンドウです。
**/
class EditAtomWindow :
    public BaseWindow
{
public:
    /* コンストラクターです。*/
    EditAtomWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
