1.
- Garbled Circuitと比べると、一ステップで加算/乗算といった論理ゲートよりも複雑な処理ができるため、同じ演算でも必要なステップ数を少なくしやすい。
- Trusted Execution Environmentと比べると、自分のハードウェア上で完結しない処理も扱うことができる。
- 加算/乗算ができるので、用途が統計処理や機械学習などに限られるDifferential PrivacyやFedearted Learningと比べると様々な計算、複雑な計算に適用できる。

2.
- アプリケーション
  - 電子投票
- 他の技術で実現可能か
  - Garbled Circuit、Secure Multi Party Computation
    - 可能（参考：https://scrapbox.io/layerx/Yao%E3%81%AEGarbled_Circuit%E3%81%AB%E3%81%A4%E3%81%84%E3%81%A6）
    - Garbled Circuitでは準同型暗号を用いる場合よりも回路が大きくなり計算コストが高くなると思われる。
    - 十分多くの人が結託した場合にセキュアでなくなる場合がある。
  - Trusted Execution Environment
    - 投票を集計するにはデータを集める必要があるが、集計サーバがTEEを使っていたところで集計サーバを信頼する必要があるため適切とは言えない。ただし、例えば集計サーバに侵入して投票データを漏洩させるような攻撃に対しては防御できると考えられる。
  - Differential Privacy
    - 各投票者の票を01とすると、これを識別できないようにするためには、集計サーバに送るデータに対し1と比べて十分大きなノイズを付加する必要がある。すべての投票を足した結果がノイズを許容できるほど大規模でかつ、数票の差を争うような投票ではないような場合には有用であると考えられる。