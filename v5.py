from typing import Dict, List
from datamodel import OrderDepth, TradingState, Order
import numpy as np


class Trader:

    def __init__(self, historical_mean=5000, window_size=17, threshold=0.01, 
                 historical_mean_spread_ratio=1, window_size_pair=100, threshold_pair=0.5, gamma=0.2):
        self.position_limits = {
            'PEARLS': 20,
            'BANANAS': 20,
            'COCONUTS': 300,
            'PINA_COLADAS': 600,
            'DIVING_GEAR': 50,
            'BERRIES': 250,
            'BAGUETTE': 150,
            'UKULELE': 70,
            'DIP': 300,
            'PICNIC_BASKET': 70

        }

        self.historical_mean = historical_mean 
        self.mean_spread_ratio = historical_mean_spread_ratio
        self.window_size = window_size 
        self.window_size_pair = window_size_pair
        self.threshold = threshold
        self.threshold_pair = threshold_pair
        self.gamma = gamma 
        self.prices_bananas = []
        self.prices_coconuts = []
        self.prices_pina_coladas = []
        self.prices_dolphin = []
        self.prices_berry = []
        self.spread_ratios = []
        self.prices_baguette = []
        self.prices_dip = []
        self.prices_ukulele = []
        self.prices_basket = []


    def run(self, state: TradingState) -> Dict[str, List[Order]]:
        """
        Only method required. It takes all buy and sell orders for all symbols as an input,
        and outputs a list of orders to be sent
        """

        # Initialize the method output dict as an empty dict
        result = {}
        alpha = 2
        LOT_SIZE = 10

        # Iterate over all the keys (the available products) contained in the order depths
        for product in state.order_depths.keys():

            # Check if the current product is the 'PEARLS' product, only then run the order logic
            if product == 'PEARLS':

                # Retrieve the Order Depth containing all the market BUY and SELL orders for PEARLS
                order_depth: OrderDepth = state.order_depths[product]

                # Initialize the list of Orders to be sent as an empty list
                orders: list[Order] = []

                # Define a fair value for the PEARLS.
                conviction = 10000 - alpha * max(state.position.get(product, 0)-10, 0) / self.position_limits[product]
                        
                # Buy orders
                if len(order_depth.sell_orders) != 0:  
                    sell_orders = order_depth.sell_orders
                    # print(sell_orders)
                    for ask, volume in sell_orders.items(): # ask volume is negative
                        if ask < conviction:
                            trade_volume = min(-volume, self.position_limits[product] - state.position.get(product, 0)) #+ve trade_volume to buy
                            
                            if trade_volume > 0:
                                orders.append(Order(product, ask, trade_volume))
                                # self.position[product] += trade_volume
                        else:
                            pass
                    
                # Sell orders
                if len(order_depth.buy_orders) != 0:
                    buy_orders = order_depth.buy_orders
                    for bid, volume in buy_orders.items(): # bid volume is positive
                        if bid > conviction:
                            trade_volume = -min(volume, self.position_limits[product] + state.position.get(product, 0)) #-ve trade_volume to sell

                            if trade_volume < 0:
                                orders.append(Order(product, bid, trade_volume))
                                # self.position[product] += trade_volume
                        else:
                            pass

                
                for order in orders:
                    if order.quantity > 0:
                        print("BUY PEAR", order.quantity, "x", order.price)
                    else:
                        print("SELL PEAR", -order.quantity, "x", order.price)

                result[product] = orders
            
            if product == 'BANANAS':
                # Retrieve the Order Depth containing all the market BUY and SELL orders for PEARLS
                order_depth: OrderDepth = state.order_depths[product]

                # Initialize the list of Orders to be sent as an empty list
                orders: list[Order] = []

                # Calculate statistics
                best_bid = max(order_depth.buy_orders.keys()) if order_depth.buy_orders else None 
                best_ask = min(order_depth.sell_orders.keys()) if order_depth.sell_orders else None

                if best_bid and best_ask:
                    mid_price = (best_bid + best_ask) // 2 
                    self.prices_bananas.append(mid_price)
                    
                
                ask_weighted = sum([p * v for p, v in order_depth.sell_orders.items()]) / sum(order_depth.sell_orders.values())
                bid_weighted = sum([p * v for p, v in order_depth.buy_orders.items()]) / sum(order_depth.buy_orders.values())
                mid_price_weighted = (ask_weighted + bid_weighted) / 2
                # price_change = self.price_lists[product][-1] - self.price_lists[product][-2] if len(self.price_lists[product]) > 1 else 0
                conviction = mid_price_weighted

                # ----
                # ACTIVE STRATEGY
                # ----

                switch = True

                # Buy orders
                if len(order_depth.sell_orders) != 0:
                    sell_orders = order_depth.sell_orders
                    # print(sell_orders)
                    for ask, volume in sell_orders.items(): # ask volume is negative
                        if ask < conviction:
                            trade_volume = min(-volume, self.position_limits[product] - state.position.get(product, 0)) #+ve trade_volume to buy
                            
                            if trade_volume > 0:
                                orders.append(Order(product, ask, trade_volume))
                                # self.position[product] += trade_volume
                                if trade_volume * 2 > self.position_limits[product] - state.position.get(product, 0):
                                    switch = False
                        else:
                            pass
                    
                # Sell orders
                if len(order_depth.buy_orders) != 0:
                    buy_orders = order_depth.buy_orders
                    for bid, volume in buy_orders.items(): # bid volume is positive
                        if bid > conviction:
                            trade_volume = -min(volume, self.position_limits[product] + state.position.get(product, 0)) #-ve trade_volume to sell

                            if trade_volume < 0:
                                orders.append(Order(product, bid, trade_volume))
                                # self.position[product] += trade_volume
                                if -trade_volume * 2 > self.position_limits[product] + state.position.get(product, 0):
                                    switch = False
                        else:
                            pass

                # ----
                # MEAN STRATEGY
                # ----

                if switch:
                    sorted_bids = sorted(order_depth.buy_orders.keys(), reverse=True)
                    sorted_asks = sorted(order_depth.sell_orders.keys())

                    # if sorted_bids and sorted_asks:
                        
                        # self.price_lists[product].append(mid_price)
                        # self.price_lists[product] = list(self.price_lists[product]).append(mid_price) 
                
                        # If there are enough data points, calculate the mean and standard deviation 
                    if len(self.prices_bananas) >= self.window_size:
                        mean = np.mean(self.prices_bananas[-self.window_size:]) 
                        stddev = np.std(self.prices_bananas[-self.window_size:])
                        position = state.position.get(product, 0)

                        for ask in sorted_asks:
                            if ask < mean - self.threshold * stddev and position < self.position_limits[product]:
                                buy_volume = min(-order_depth.sell_orders[ask], self.position_limits[product] - position)
                                orders.append(Order(product, ask, buy_volume))
                                # self.position[product] += buy_volume
                        
                        for bid in sorted_bids:
                            if bid > mean + self.threshold * stddev and position > -self.position_limits[product]:
                                sell_volume = min(order_depth.buy_orders[bid], self.position_limits[product] + position)
                                orders.append(Order(product, bid, -sell_volume))
                                # self.position[product] -= sell_volume

                for order in orders:
                    if order.quantity > 0:
                        print("BUY",product, order.quantity, "x", order.price)
                    else:
                        print("SELL",product, -order.quantity, "x", order.price)

                result[product] = orders
            
            if product == 'COCONUTS':
                order_depth_coconuts: OrderDepth = state.order_depths[product]
                order_depth_pina: OrderDepth = state.order_depths['PINA_COLADAS']

                orders_coconuts: list[Order] = []

                best_bid_coco = max(order_depth_coconuts.buy_orders.keys()) if order_depth_coconuts.buy_orders.keys() else None 
                best_ask_coco = min(order_depth_coconuts.sell_orders.keys()) if order_depth_coconuts.sell_orders.keys() else None
                if best_bid_coco and best_ask_coco:
                    mid_price_coco = (best_bid_coco + best_ask_coco) // 2 

                best_bid_pina = max(order_depth_pina.buy_orders.keys()) if order_depth_pina.buy_orders.keys() else None 
                best_ask_pina = min(order_depth_pina.sell_orders.keys()) if order_depth_pina.sell_orders.keys() else None
                if best_bid_pina and best_ask_pina:
                    mid_price_pina = (best_bid_pina + best_ask_pina) // 2 

                self.prices_coconuts.append(mid_price_coco)
                self.prices_pina_coladas.append(mid_price_pina)
                spread_ratio = mid_price_pina / mid_price_coco
                self.spread_ratios.append(spread_ratio)
                self.mean_spread_ratio = self.gamma * spread_ratio + (1 - self.gamma) * self.mean_spread_ratio

                # If there are enough data points, calculate the standard deviation of the spread ratios
                if len(self.spread_ratios) >= self.window_size_pair:
                    stddev = np.std(self.spread_ratios[-self.window_size_pair:])

                    # If the spread ratio is below the lower threshold, consider buying coconuts and selling pina coladas
                    if spread_ratio > self.mean_spread_ratio - self.threshold_pair * stddev:
                        remaining_coconuts = self.position_limits[product] - state.position.get(product, 0) #+ve

                        if remaining_coconuts > 0:
                            buy_volume = min(remaining_coconuts, LOT_SIZE) # >0
                            ask = list(order_depth_coconuts.sell_orders.keys())[0] -1
                            orders_coconuts.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume

                    # If the spread ratio is above the upper threshold, consider selling coconuts and buying pina coladas
                    elif spread_ratio < self.mean_spread_ratio + self.threshold_pair * stddev:
                        remaining_coconuts = state.position.get(product, 0)
                        
                        if remaining_coconuts > 0:
                            sell_volume = min(remaining_coconuts, LOT_SIZE)
                            bid = list(order_depth_coconuts.buy_orders.keys())[-1] +1
                            orders_coconuts.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= sell_volume

                for order in orders_coconuts:
                    if order.quantity > 0:
                        print("BUY COCO", order.quantity, "x", order.price)
                    else:
                        print("SELL COCO", -order.quantity, "x", order.price)

                result[product] = orders_coconuts
            
            """
            if product == 'PINA_COLADAS':
                order_depth_coconuts: OrderDepth = state.order_depths['COCONUTS']
                order_depth_pina: OrderDepth = state.order_depths['PINA_COLADAS']

                orders_pina_coladas: list[Order] = []

                best_bid_coco = max(order_depth_coconuts.buy_orders.keys()) if order_depth_coconuts.buy_orders.keys() else None 
                best_ask_coco = min(order_depth_coconuts.sell_orders.keys()) if order_depth_coconuts.sell_orders.keys() else None
                if best_bid_coco and best_ask_coco:
                    mid_price_coco = (best_bid_coco + best_ask_coco) // 2 

                best_bid_pina = max(order_depth_pina.buy_orders.keys()) if order_depth_pina.buy_orders.keys() else None 
                best_ask_pina = min(order_depth_pina.sell_orders.keys()) if order_depth_pina.sell_orders.keys() else None
                if best_bid_pina and best_ask_pina:
                    mid_price_pina = (best_bid_pina + best_ask_pina) // 2 

                self.prices_coconuts.append(mid_price_coco)
                self.prices_pina_coladas.append(mid_price_pina)
                spread_ratio = mid_price_pina / mid_price_coco
                self.spread_ratios.append(spread_ratio)
                self.mean_spread_ratio = self.gamma * spread_ratio + (1 - self.gamma) * self.mean_spread_ratio

                # If there are enough data points, calculate the standard deviation of the spread ratios
                if len(self.spread_ratios) >= self.window_size_pair:
                    stddev = np.std(self.spread_ratios[-self.window_size_pair:])

                    # If the spread ratio is below the lower threshold, consider buying coconuts and selling pina coladas
                    if spread_ratio < self.mean_spread_ratio - self.threshold_pair * stddev:
                        remaining_pina_coladas = state.position.get(product, 0)

                        if remaining_pina_coladas > 0:
                            sell_volume = min(remaining_pina_coladas, 1) #>0
                            bid = list(order_depth_pina.buy_orders.keys())[-1] + 1
                            orders_pina_coladas.append(Order(product, bid, -sell_volume))
                            state[product] -= sell_volume

                    # If the spread ratio is above the upper threshold, consider selling coconuts and buying pina coladas
                    elif spread_ratio > self.mean_spread_ratio + self.threshold_pair * stddev:
                        remaining_pina_coladas = self.position_limits[product] - state.position.get(product, 0)
              
                        if remaining_pina_coladas > 0:
                            buy_volume = min(remaining_pina_coladas, 1)
                            ask = list(order_depth_pina.sell_orders.keys())[0]
                            orders_pina_coladas.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume

                for order in orders_pina_coladas:
                        if order.quantity > 0:
                            print("BUY PINA", order.quantity, "x", order.price)
                        else:
                            print("SELL PINA", -order.quantity, "x", order.price)

                result[product] = orders_pina_coladas
            """

            if product == 'DIVING_GEAR':
                order_depth: OrderDepth = state.order_depths[product]
                orders: list[Order] = []
                self.prices_dolphin.append(state.observations['DOLPHIN_SIGHTINGS'])

                if len(self.prices_dolphin) < 30:
                    ask = list(order_depth.sell_orders.keys())[0] - 1
                    buy_volume = min(LOT_SIZE, self.position_limits[product] - state.position.get(product, 0))
                    orders.append(Order(product, ask, buy_volume))


                elif len(self.prices_dolphin) >= 30:
                    short_mean = np.mean(self.prices_dolphin[-5:]) 
                    long_mean = np.mean(self.prices_dolphin[-30:]) 
                
                    if short_mean > 3 * long_mean:
                        if state.position.get(product, 0) < self.position_limits[product]:
                            ask = list(order_depth.sell_orders.keys())[0] - 1
                            buy_volume = min(LOT_SIZE, self.position_limits[product] - state.position.get(product, 0))
                            orders.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume
                        
                    elif short_mean < 3 * long_mean:
                        if state.position.get(product, 0) > -self.position_limits[product]:
                            bid = list(order_depth.buy_orders.keys())[-1] + 1
                            sell_volume = min(LOT_SIZE, self.position_limits[product] + state.position.get(product, 0)) 
                            orders.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= sell_volume

                for order in orders:
                        if order.quantity > 0:
                            print("BUY DIVE", order.quantity, "x", order.price)
                        else:
                            print("SELL DIVE", -order.quantity, "x", order.price)
                result[product] = orders

            """
            if product == 'BERRIES':
                order_depth: OrderDepth = state.order_depths[product]
                orders: list[Order] = []

                time: int = state.timestamp
                mid_day = 5000
                buffer = 100*20

                # Before and after mid-day, do market making
                if time < mid_day - buffer or time > mid_day + buffer:
                    # Calculate statistics
                    best_bid = max(order_depth.buy_orders.keys()) if order_depth.buy_orders else None 
                    best_ask = min(order_depth.sell_orders.keys()) if order_depth.sell_orders else None

                    if best_bid and best_ask:
                        mid_price = (best_bid + best_ask) // 2 
                        self.prices_berry.append(mid_price)
                    
                    ask_weighted = sum([p * v for p, v in order_depth.sell_orders.items()]) / sum(order_depth.sell_orders.values())
                    bid_weighted = sum([p * v for p, v in order_depth.buy_orders.items()]) / sum(order_depth.buy_orders.values())
                    mid_price_weighted = (ask_weighted + bid_weighted) / 2
                    conviction = mid_price_weighted

                    # ----
                    # ACTIVE STRATEGY
                    # ----

                    switch = True

                    # Buy orders
                    if len(order_depth.sell_orders) != 0:
                        sell_orders = order_depth.sell_orders
                        # print(sell_orders)
                        ask = list(sell_orders.keys())[0]-1 # ask volume is negative
                        if ask < conviction:
                            trade_volume = min(LOT_SIZE, self.position_limits[product] - state.position.get(product, 0)) #+ve trade_volume to buy
                            if trade_volume > 0:
                                orders.append(Order(product, ask, trade_volume))
                                if trade_volume * 2 > self.position_limits[product] - state.position.get(product, 0):
                                    switch = False
                        else:
                            pass
                        
                    # Sell orders
                    if len(order_depth.buy_orders) != 0:
                        buy_orders = order_depth.buy_orders
                        bid = list(buy_orders.keys())[-1]+1 # bid volume is positive
                        if bid > conviction:
                            trade_volume = -min(LOT_SIZE, self.position_limits[product] + state.position.get(product, 0)) #-ve trade_volume to sell
                            if trade_volume < 0:
                                orders.append(Order(product, bid, trade_volume))
                                if -trade_volume * 2 > self.position_limits[product] + state.position.get(product, 0):
                                    switch = False
                        else:
                            pass

                    # ----
                    # MEAN STRATEGY
                    # ----

                    if switch:
                        sorted_bids = sorted(order_depth.buy_orders.keys(), reverse=True)
                        sorted_asks = sorted(order_depth.sell_orders.keys())

                        if sorted_bids and sorted_asks:
                            
                            self.prices_berry.append(mid_price)
                    
                            # If there are enough data points, calculate the mean and standard deviation 
                            if len(self.prices_berry) >= self.window_size:
                                mean = np.mean(self.prices_berry[-self.window_size:]) 
                                stddev = np.std(self.prices_berry[-self.window_size:])
                                position = state.position.get(product, 0)

                                for ask in sorted_asks:
                                    if ask < mean - self.threshold * stddev and position < self.position_limits[product]:
                                        buy_volume = min(-order_depth.sell_orders[ask], self.position_limits[product] - position)
                                        orders.append(Order(product, ask, buy_volume))
                                
                                for bid in sorted_bids:
                                    if bid > mean + self.threshold * stddev and position > -self.position_limits[product]:
                                        sell_volume = min(order_depth.buy_orders[bid], self.position_limits[product] + position)
                                        orders.append(Order(product, bid, -sell_volume))
                    
                # Close to mid-day, fill up position
                elif time > mid_day - buffer and time <= mid_day:
                    position = state.position.get(product, 0)
                    
                    # As long as there is a position, accept best ask offer
                    if position < self.position_limits[product]:
                        best_ask = min(order_depth.sell_orders.keys()) if order_depth.sell_orders else None
                        if best_ask:
                            buy_volume = min(-order_depth.sell_orders[best_ask], self.position_limits[product] - position)
                            orders.append(Order(product, best_ask, buy_volume))

                # After mid-day but close to it, fill up negative position
                elif time > mid_day and time < mid_day + buffer:
                    position = state.position.get(product, 0)
                    
                    # As long as there is a position, accept best bid offer
                    if position > -self.position_limits[product]:
                        best_bid = max(order_depth.buy_orders.keys()) if order_depth.buy_orders else None
                        if best_bid:
                            sell_volume = min(order_depth.buy_orders[best_bid], self.position_limits[product] + position)
                            orders.append(Order(product, best_bid, -sell_volume))

                else:
                    pass

                for order in orders:
                        if order.quantity > 0:
                            print("BUY BERY", order.quantity, "x", order.price)
                        else:
                            print("SELL BERY", -order.quantity, "x", order.price)

                result[product] = orders
            """
            
            if product == 'BAGUETTE':
                # Retrieve the Order Depth containing all the market BUY and SELL orders for PEARLS
                order_depth: OrderDepth = state.order_depths[product]

                # Initialize the list of Orders to be sent as an empty list
                orders: list[Order] = []

                ask_weighted = sum([p * v for p, v in order_depth.sell_orders.items()]) / sum(order_depth.sell_orders.values())
                bid_weighted = sum([p * v for p, v in order_depth.buy_orders.items()]) / sum(order_depth.buy_orders.values())
                mid_price_weighted = (ask_weighted + bid_weighted) // 2
                self.prices_baguette.append(mid_price_weighted)

                if len(self.prices_baguette) >= 50:
                    short_mean = np.mean(self.prices_baguette[-5:]) 
                    long_mean = np.mean(self.prices_baguette[-50:]) 
                
                    if short_mean > long_mean + 0.5:
                        if state.position.get(product, 0) < self.position_limits[product]:
                            ask = list(order_depth.sell_orders.keys())[0] - 2
                            buy_volume = min(5, self.position_limits[product] - state.position.get(product, 0))
                            orders.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume
                        
                    elif short_mean < long_mean - 0.5:
                        if state.position.get(product, 0) > -self.position_limits[product]:
                            bid = list(order_depth.buy_orders.keys())[-1] + 2
                            sell_volume = min(5, self.position_limits[product] + state.position.get(product, 0)) 
                            orders.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= sell_volume

                for order in orders:
                    if order.quantity > 0:
                        print("BUY",product, order.quantity, "x", order.price)
                    else:
                        print("SELL",product, -order.quantity, "x", order.price)

                result[product] = orders

            if product == 'DIP':
                # Retrieve the Order Depth containing all the market BUY and SELL orders for PEARLS
                order_depth: OrderDepth = state.order_depths[product]

                # Initialize the list of Orders to be sent as an empty list
                orders: list[Order] = []

                ask_weighted = sum([p * v for p, v in order_depth.sell_orders.items()]) / sum(order_depth.sell_orders.values())
                bid_weighted = sum([p * v for p, v in order_depth.buy_orders.items()]) / sum(order_depth.buy_orders.values())
                mid_price_weighted = (ask_weighted + bid_weighted) // 2
                self.prices_dip.append(mid_price_weighted)

                if len(self.prices_dip) >= 40:
                    short_mean = np.mean(self.prices_dip[-5:]) 
                    long_mean = np.mean(self.prices_dip[-40:]) 
                
                    if short_mean > long_mean+0.2:
                        if state.position.get(product, 0) < self.position_limits[product]:
                            ask = list(order_depth.sell_orders.keys())[0] - 1
                            buy_volume = min(LOT_SIZE, self.position_limits[product] - state.position.get(product, 0))
                            orders.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume
                        
                    elif short_mean < long_mean-0.2:
                        if state.position.get(product, 0) > -self.position_limits[product]:
                            bid = list(order_depth.buy_orders.keys())[-1] + 1
                            sell_volume = min(LOT_SIZE, self.position_limits[product] + state.position.get(product, 0)) 
                            orders.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= trade_volume

                for order in orders:
                    if order.quantity > 0:
                        print("BUY",product, order.quantity, "x", order.price)
                    else:
                        print("SELL",product, -order.quantity, "x", order.price)

                result[product] = orders

            if product == 'UKULELE':
                # Retrieve the Order Depth containing all the market BUY and SELL orders for PEARLS
                order_depth: OrderDepth = state.order_depths[product]

                # Initialize the list of Orders to be sent as an empty list
                orders: list[Order] = []

                ask_weighted = sum([p * v for p, v in order_depth.sell_orders.items()]) / sum(order_depth.sell_orders.values())
                bid_weighted = sum([p * v for p, v in order_depth.buy_orders.items()]) / sum(order_depth.buy_orders.values())
                mid_price_weighted = (ask_weighted + bid_weighted) // 2
                self.prices_ukulele.append(mid_price_weighted)

                if len(self.prices_ukulele) >= 40:
                    short_mean = np.mean(self.prices_ukulele[-5:]) 
                    long_mean = np.mean(self.prices_ukulele[-40:]) 
                
                    if short_mean > long_mean+0.2:
                        if state.position.get(product, 0) < self.position_limits[product]:
                            ask = list(order_depth.sell_orders.keys())[0] - 2
                            buy_volume = min(5, self.position_limits[product] - state.position.get(product, 0))
                            orders.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume
                        
                    elif short_mean < long_mean-0.2:
                        if state.position.get(product, 0) > -self.position_limits[product]:
                            bid = list(order_depth.buy_orders.keys())[-1] + 2
                            sell_volume = min(5, self.position_limits[product] + state.position.get(product, 0)) 
                            orders.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= sell_volume

                for order in orders:
                    if order.quantity > 0:
                        print("BUY",product, order.quantity, "x", order.price)
                    else:
                        print("SELL",product, -order.quantity, "x", order.price)

                result[product] = orders

            if product == 'PICNIC_BASKET':
                if len(self.prices_baguette) >= 1 and len(self.prices_dip) >= 1 and len(self.prices_ukulele) >= 1:
                    mid_price_weighted = (2*self.prices_baguette[-1] + self.prices_dip[-1] + self.prices_ukulele[-1]) // 4
                    self.prices_basket.append(mid_price_weighted)

                if len(self.prices_basket) >= 40:
                    short_mean = np.mean(self.prices_basket[-5:]) 
                    long_mean = np.mean(self.prices_basket[-40:]) 
                
                    if short_mean > long_mean+0.05:
                        if state.position.get(product, 0) < self.position_limits[product]:
                            ask = list(order_depth.sell_orders.keys())[0]
                            buy_volume = min(LOT_SIZE, self.position_limits[product] - state.position.get(product, 0))
                            orders.append(Order(product, ask, buy_volume))
                            # self.position[product] += buy_volume
                        
                    elif short_mean < long_mean-0.05:
                        if state.position.get(product, 0) > -self.position_limits[product]:
                            bid = list(order_depth.buy_orders.keys())[-1]
                            sell_volume = min(LOT_SIZE, self.position_limits[product] + state.position.get(product, 0)) 
                            orders.append(Order(product, bid, -sell_volume))
                            # self.position[product] -= sell_volume

                for order in orders:
                    if order.quantity > 0:
                        print("BUY",product, order.quantity, "x", order.price)
                    else:
                        print("SELL",product, -order.quantity, "x", order.price)


        # Print position
        try:
            print("POSITION", state.position)
        except:
            print("POSITION", 0)

        return result